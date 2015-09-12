from __future__ import print_function

import os
import time
import sqlite3
import gzip
import threading
from multiprocessing import Queue, Process
import logging
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s|%(levelname)s|%(message)s')

import numpy as np
import pandas as pd

from utils import pretty_usage, timeit, memo, split
from objs import BloomFilterBuilder


PROB = 1e-9
K_MER_SIZE = 20
NUM_CPUS = 32         # given 32 cores

def generate_bf_seqs():
    """generate subarray of sequences needed to create a bloomfilter"""
    input_ = '/projects/btl2/zxue/microorganism_profiling/libraries_assesment/gg_unite/gg_unite.csv.gz'
    logging.info('reading {0}'.format(input_))
    df = pd.read_csv(input_, compression='gzip')
    # df = pd.read_csv(input_, compression='gzip', nrows=50)
    logging.info('reading Done')

    seqs = df.seq.values
    logging.info('Given {0} sequences: about 2x bloomfilters will be '
                 'generated'.format(df.seq.shape[0]))

    # index each seq, 0-based index, ends up with tuples of (seq_id, seq)
    id_seqs = zip(xrange(len(seqs)), seqs)

    for ele in split(id_seqs):
        # each element is of structure:
        # (bf_id, [(seq_id1, seq1), (seq_id2, seq2), ..])
        yield ele

@memo
def kmerize(seq, kmer_size):
    res = []
    for i in xrange(len(seq) + 1 - kmer_size):
        res.append(seq[i:i + kmer_size])
    return res


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


def calc_fp(m ,h, n):
    """
    Calculate false positive rate

    :param m: size of bloomfilter
    :param h: number of hash count
    :param n: number of elements

    ref: http://www.maxburstein.com/blog/creating-a-simple-bloom-filter/
    """
    return (1 - (1 - 1 / m) ** (h * n)) ** h

# @timeit
def generate_bf(id_seqs):
    if len(id_seqs) == 1:
        seq_id = id_seqs[0][0]
    else:
        seq_id = None

    num_total_kmers = 0                 # include all kmers
    uniq_kmers = []
    # _: don't track id if there are multiple sequences
    for (_, seq) in id_seqs:
        kmers = kmerize(seq, K_MER_SIZE)
        num_total_kmers += len(kmers)
        # intermediate set to reduce memory usage
        uniq_kmers.extend(list(set(kmers)))

    uniq_kmers = set(uniq_kmers)    # final set
    num_uniq_kmers = len(uniq_kmers)

    # logging.info(num_total_kmers, num_uniq_kmers)

    hash_count, bf_size = calc_hash_count_and_bf_size(num_uniq_kmers, PROB)
    # logging.info(hash_count)
    # logging.info(bf_size, pretty_usage(bf_size / 8.))
    # logging.info("real false positive rate: {0}".format(calc_fp(bf_size, hash_count, num_uniq_kmers)))

    bf_builder = BloomFilterBuilder(bf_size, hash_count)
    for kmer in uniq_kmers:
        bf_builder.add(kmer)

    return bf_builder.bit_array, seq_id


def worker(queue, db_count):
    logging.info('connecting to db...')
    output_db = "db/{0}.db".format(db_count)
    if output_db:
        os.remove(output_db)
    conn = sqlite3.connect(output_db)
    cursor = conn.cursor()

    # cursor.execute("CREATE TABLE bloomfilter (bf_id INTEGER PRIMARY KEY, bitarray BLOB, FOREIGN KEY(seq_id) REFERENCES seq(id))")
    cursor.execute("CREATE TABLE bloomfilter (bf_id INTEGER PRIMARY KEY, bitarray BLOB, seq_id)")

    counter = 0
    while True:
        bf_id, id_seqs = queue.get()
        bf, seq_id = generate_bf(id_seqs)
        if seq_id:
            cursor.execute("INSERT INTO bloomfilter values (?, ?, ?)",
                           (bf_id, sqlite3.Binary(bf.tobytes()), seq_id))
        else:
            cursor.execute("INSERT INTO bloomfilter values (?, ?, ?)",
                           (bf_id, sqlite3.Binary(bf.tobytes()), None))
        counter += 1
        if counter % 500 == 0:
            logging.info("commit {0} executions".format(counter))
            conn.commit()

        if queue.empty():
            logging.info('db_count: {0}: waiting 5s for queue to be non-empty'.format(db_count))
            time.sleep(5)
            if queue.empty():   # still empty
                logging.info('db_count: {0}: waiting 5s for queue to be non-empty'.format(db_count))
                break
    conn.commit()               # commit remaining execution


@timeit
def put_to_db(res_queue, output):
    logging.info('connecting to db...')
    conn = sqlite3.connect(":memory:")
    cursor = conn.cursor()

    # cursor.execute("CREATE TABLE bloomfilter (bf_id INTEGER PRIMARY KEY, bitarray BLOB, FOREIGN KEY(seq_id) REFERENCES seq(id))")
    cursor.execute("CREATE TABLE bloomfilter (bf_id INTEGER PRIMARY KEY, bitarray BLOB, seq_id)")

    while not res_queue.empty():
        bf_id, bf, seq_id = res_queue.get()
        # logging.info('putting {0}th bf with seq_id {1} into db'.format(bf_id, seq_id))
        if seq_id:
            cursor.execute("INSERT INTO bloomfilter values (?, ?, ?)",
                           (bf_id, sqlite3.Binary(bf.tobytes()), seq_id))
        else:
            cursor.execute("INSERT INTO bloomfilter values (?, ?, ?)",
                           (bf_id, sqlite3.Binary(bf.tobytes()), None))
        res_queue.task_done()

    logging.info('committing db')
    conn.commit()
    logging.info('dumping db to {0}'.format(output))
    dump_db(conn, output)
    conn.close()


@timeit
def dump_db(conn, output):
    with gzip.open(output, 'w') as opf:
        for line in conn.iterdump():
            opf.write('%s\n' % line)

    
def main():
    iter_bf_seqs = generate_bf_seqs()

    queue = Queue()

    procs = []
    for i in range(NUM_CPUS - 1):
        proc = Process(target=worker, args=(queue, i))
        proc.daemon = True
        procs.append(proc)
        proc.start()

    for k, (bf_id, id_seqs) in enumerate(iter_bf_seqs):
        if (k + 1) % 10000 == 0:
            logging.info('enqueuing {0}th bf'.format(k + 1))
        queue.put((bf_id, id_seqs))
    logging.info('enqueued {0}th bf'.format(k + 1))

    for proc in procs:
        proc.join()

    # output = 'dump.sql'
    # put_to_db(res_queue, output)



if __name__ == "__main__":
    main()
