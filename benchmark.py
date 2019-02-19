from __future__ import print_function

import sqlite3
import gzip
import threading
import Queue
import logging
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s|%(levelname)s|%(message)s')

import numpy as np
import pandas as pd

from utils import pretty_usage, timeit, split
from objs import BloomFilterBuilder


PROB = 1e-9
K_MER_SIZE = 20
NUM_WORKER_THREADS = 32         # given 32 cores


def generate_bf_seqs():
    """generate subarray of sequences needed to create a bloomfilter"""
    df = pd.read_csv('/projects/btl2/zxue/microorganism_profiling/libraries_assesment/gg_unite/gg_unite.csv.gz',
                     # compression='gzip', nrows=50)
                     compression='gzip')

    seqs = df.seq.values
    logging.info('Given {0} sequences: about 2x bloomfilters will be '
                 'generated'.format(df.seq.shape[0]))

    # index each seq, 0-based index, ends up with tuples of (seq_id, seq)
    id_seqs = zip(xrange(len(seqs)), seqs)

    for ele in split(id_seqs):
        # each element is of structure:
        # (bf_id, [(seq_id1, seq1), (seq_id2, seq2), ..])
        yield ele


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


@timeit
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


def worker(queue, res_queue):
    while True:
        bf_id, id_seqs = queue.get()
        bf, seq_id = generate_bf(id_seqs)
        res_queue.put((bf_id, bf, seq_id))
        queue.task_done()


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

    
def main(output):
    queue = Queue.Queue()
    res_queue = Queue.Queue()   # Queue for puting results

    for i in range(NUM_WORKER_THREADS - 1):
        thr = threading.Thread(target=worker, args=(queue, res_queue))
        thr.daemon = True
        thr.start()

    # for (bf_id, id_seqs) in generate_bf_seqs():
    #     queue.put((bf_id, id_seqs))
    for k, (bf_id, id_seqs) in enumerate(generate_bf_seqs()):
        if k % 50 == 0:
            logging.info('enqueuing {0}th bf'.format(k))
        queue.put((bf_id, id_seqs))

    queue.join()                # block till all tasks are finished

    put_to_db(res_queue, output)
    res_queue.join()


if __name__ == "__main__":
    # output = 'lele.sql'
    # main(output)
    for k, (bf_id, id_seqs) in enumerate(generate_bf_seqs()):
        print(bf_id)
        # if k % 50 == 0:
        #     logging.info('enqueuing {0}th bf'.format(k))
