
from __future__ import print_function

import glob
import os
import time
import sqlite3
import gzip
import threading
from multiprocessing import Queue, Process
import logging

import numpy as np
import pandas as pd

from utils import pretty_usage, timeit, memo, split_indexes, kmerize, calc_level
from objs import BloomFilterBuilder

from settings import DEBUG


from constants import (
    K_MER_SIZE,
    NBR)

if DEBUG:
    NUM_CPUS = 4
    DB_OUTPUT_DIR = 'debug_db_nbr_{0}'.format(NBR)
else:
    NUM_CPUS = 32         # given 32 cores
    DB_OUTPUT_DIR = 'db_nbr_{0}'.format(NBR)
if not os.path.exists(DB_OUTPUT_DIR):
    os.mkdir(DB_OUTPUT_DIR)



if DEBUG:
    logging.basicConfig(level=logging.DEBUG,
                        format='%(asctime)s|%(levelname)s|%(message)s')
else:
    logging.basicConfig(level=logging.DEBUG,
                        filename='build_tree_nbr_{0}.log'.format(NBR), filemode='w',
                        format='%(asctime)s|%(levelname)s|%(message)s')


SQL_CREATE_TABLE = """CREATE TABLE 
    bloomfilter
    (bf_id INTEGER PRIMARY KEY,
     level INTEGER,
    num_uniq_kmers INTEGER,
    size INTEGER,
    hash_count INTEGER,
    fpr REAL,
    bitarray BLOB,
    seq_id)"""


def generate_seqid_seqs():
    """return "a list of tuples of (seqid, seq)"""

    input_ = '/projects/btl2/zxue/microorganism_profiling/libraries_assesment/gg_unite/gg_unite.csv.gz'
    logging.info('reading {0}'.format(input_))
    if DEBUG:
        df = pd.read_csv(input_, compression='gzip', nrows=1000)
    else:
        df = pd.read_csv(input_, compression='gzip')
    logging.info('reading Done')

    seqs = df.seq.values
    logging.info('Given {0} sequences: about 2x bloomfilters will be '
                 'generated'.format(df.seq.shape[0]))

    # index each seq, 0-based index, ends up with tuples of (seq_id, seq)
    seqid_seqs = zip(xrange(len(seqs)), seqs)
    return seqid_seqs


def calc_bf_size(n, p):
    return - (n * np.log(p)) / (np.log(2) ** 2)


def calc_hash_count_and_bf_size(n, p):
    """
    :param n: number of elements
    :param p: false positive probability

    ref: http://www.maxburstein.com/blog/creating-a-simple-bloom-filter/
    """
    h = - np.log(p) / np.log(2)  # hash_count
    m = n * h / np.log(2)      # bloomfilter size
    return int(h), int(m)


def build_bf(id_seqs):
    """just estimate number of uniq k_kmers to speed up processing"""
    if len(id_seqs) == 1:
        seq_id = id_seqs[0][0]
    else:
        seq_id = None

    calc_num_kmers = lambda seq: (len(seq) - K_MER_SIZE + 1)
    num_kmers = sum(calc_num_kmers(x[1]) for x in id_seqs)

    prob = 0.0075

    hash_count, bf_size = calc_hash_count_and_bf_size(num_kmers, prob)
    if DEBUG:
        logging.info('bf_size: {0} ({1})'.format(bf_size, pretty_usage(bf_size / 8.)))
        logging.info('hash_count: {0}'.format(hash_count))

    # The number of unique kmers is approximate because it didn't take into
    # account of same kmers amond different sequences. This is for saving
    # memory and speed
    num_total_kmers = 0
    bf_builder = BloomFilterBuilder(bf_size, hash_count)
    bf_add = bf_builder.add_uniq
    # _: don't track id if there are multiple sequences
    for k, (_, seq) in enumerate(id_seqs):
        if (k + 1) % 10000 == 0:
            logging.info('working on {0}th seqid_seqs'.format(k + 1))
        kmers = kmerize(seq, K_MER_SIZE)
        num_total_kmers += len(kmers)
        kmers_set = set(kmers)
        for __ in kmers_set:
            bf_add(__)

    num_uniq_kmers = bf_builder.num_uniq_elems
    fpr = bf_builder.calc_fpr()

    if DEBUG:
        logging.info('Total number of kmers: {0}'.format(num_total_kmers))
        logging.info("Number of uniq kmers: {0}".format(num_uniq_kmers))
        logging.info("Real false positive rate: {0}".format(fpr))

    return (num_uniq_kmers, bf_builder.size, bf_builder.hash_count,
            fpr, bf_builder.bit_array, seq_id)


def worker(queue, db_count):
    logging.info('connecting to db...')
    output_db = os.path.join(DB_OUTPUT_DIR, "{0}.db".format(db_count))
    if os.path.exists(output_db):
        os.remove(output_db)
    conn = sqlite3.connect(output_db)
    cursor = conn.cursor()

    cursor.execute(SQL_CREATE_TABLE)

    counter = 0
    while True:
        bf_id, id_seqs = queue.get()
        level = calc_level(bf_id, NBR)
        num_uniq_kmers, bf_size, bf_hash_count, bf_fpr, bf, seq_id = build_bf(id_seqs)
        # seq_id is None if the bf is built on multiple sequences
        cursor.execute("INSERT INTO bloomfilter values (?, ?, ?, ?, ?, ?, ?, ?)",
                       (bf_id, level, num_uniq_kmers, bf_size, bf_hash_count, bf_fpr,
                        sqlite3.Binary(bf.tobytes()), seq_id))
        counter += 1
        if counter % 10000 == 0:
            logging.info("commit {0} executions".format(counter))
            conn.commit()

        if queue.empty():
            more_task = False   # indicate whether there are more tasks or not
            if DEBUG:
                n_loop = 3
            else:
                n_loop = 10
            for i in xrange(1, n_loop + 1):
                logging.info('db_count: {0}: waiting {1}s for queue to be non-empty'.format(db_count, i))
                time.sleep(i)
                if not queue.empty():
                    logging.info('db_count: {0}: more tasks found after {1}s '
                                 'wait, back to work'.format(db_count, (i + 1) * i / 2))
                    more_task = True
                    break

            if more_task:
                continue
            else:
                total_seconds = (n_loop + 1) * n_loop / 2.
                logging.info('db_count: {0}: still empty after {1}s wait for '
                             'the queue, break the while loop'.format(
                                 db_count, total_seconds))
                break
    conn.commit()               # commit remaining execution



def combine_db(db_dir):
    output_db = os.path.join(db_dir, 'combined.db')

    if os.path.exists(output_db):
        os.remove(output_db)


    dbs = sorted(glob.glob(os.path.join(db_dir, '[0-9]*.db')))

    conn = sqlite3.connect(output_db)
    cursor = conn.cursor()

    cursor.execute(SQL_CREATE_TABLE)

    for db in dbs:
        logging.info('combining {0} into {1}'.format(db, output_db))
        cursor.execute('attach "{0}" as toCombine'.format(db))
        cursor.execute('insert into bloomfilter select * from toCombine.bloomfilter'.format(db))
        cursor.execute('detach toCombine'.format(db))

    conn.commit()
    conn.close()


def main():
    seqid_seqs = generate_seqid_seqs()
    bfid_beg_end_generator = split_indexes(0, len(seqid_seqs), nbr=NBR)

    queue = Queue()

    procs = []
    for i in range(NUM_CPUS - 1):
        proc = Process(target=worker, args=(queue, i))
        proc.daemon = True
        procs.append(proc)
        proc.start()

    freq, next_freq = 1, 10
    for k, (bf_id, beg, end) in enumerate(bfid_beg_end_generator):
        k_plus_1 = k + 1
        if k_plus_1 <= next_freq:
            if k_plus_1 % freq == 0:
                logging.info('enqueuing {0}th bf'.format(k + 1))
        else:
            freq *= 10
            next_freq = freq * 10

        queue.put((bf_id, seqid_seqs[beg: end]))
    logging.info('enqueued {0}th bf'.format(k + 1))

    for proc in procs:
        proc.join()

    combine_db(DB_OUTPUT_DIR)


if __name__ == "__main__":
    main()
