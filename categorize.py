from __future__ import unicode_literals

import os
import time
import bitarray
import sqlite3
import pandas as pd
from multiprocessing import Manager, Value, Queue, Process, Lock
import gzip

from Bio import SeqIO

from objs import BloomFilter
import utils as U
from constants import (
    K_MER_SIZE,
    SCORE_CUTOFF,
    NBR,
    COMPLEMENT_DD)

import logging

from settings import DEBUG

if DEBUG:
    logging.basicConfig(level=logging.DEBUG,
                        format='%(asctime)s|%(levelname)s|%(message)s')
else:
    logging.basicConfig(level=logging.DEBUG,
                        filename='categorize.log', filemode='w',
                        format='%(asctime)s|%(levelname)s|%(message)s')

if DEBUG:
    NUM_CPUS = 4
    DB_OUTPUT_DIR = 'debug_db'
else:
    NUM_CPUS = 4         # given 32 cores
    DB_OUTPUT_DIR = 'db'
if not os.path.exists(DB_OUTPUT_DIR):
    os.mkdir(DB_OUTPUT_DIR)


NO_MATCH = 'no_match.txt'


def load_bfs(db_file):
    """
    :param db_file: a sqlite3 database file all bfs stored
    """

    logging.info('loading {0} into memory'.format(db_file))
    conn = sqlite3.connect(db_file)
    cur = conn.cursor()
    
    sql = 'select bf_id, size, hash_count, bitarray from bloomfilter'
    logging.info('executing {0}'.format(sql))
    cur.execute(sql)
    
    logging.info('building bfs')
    bfs = {}
    freq = 1
    next_freq = freq * 10
    for k, (bf_id, size, hash_count, bf) in enumerate(cur):
        bfs[bf_id] = BloomFilter(size, hash_count, bf)
        k_plus_1 = k + 1
        if k_plus_1 <= next_freq:
            if k_plus_1 % freq == 0:
                logging.info('built {0} bfs'.format(k_plus_1))
        else:
            freq *= 10
            next_freq = freq * 10
    logging.info('built {0}th bfs in total'.format(k_plus_1))
    conn.close()

    logging.info('loading is done'.format(db_file))
    return bfs


def reverse_complement(seq):
    return ''.join([COMPLEMENT_DD.get(x, 'N') for x in seq[::-1]])


def calc_score(read, k_mer_size, bf):
    l = len(read)
    k = k_mer_size
    d = float(l - k)

    s = 0
    j = 0
    for kmer in U.kmerize(read, k):
        if kmer in bf:
            j += 1
            s += (1 - 1 / (j + 1)) / d
        else:
            j = 0
    return s


def score(read, bf):
    f_read = read
    r_read = reverse_complement(f_read)

    f_s = calc_score(f_read, K_MER_SIZE, bf)
    r_s = calc_score(r_read, K_MER_SIZE, bf)

    # print f_s, r_s
    return max(f_s, r_s)


def get_bf(read, bfs, nbr, score_cutoff, hit, bf_id=0, level=0):
    """
    if calculated score < score_cutoff, then stop searching

    :param hit: a list that whole hit bf_ids
    """

    cids = U.calc_children_id(bf_id, level, nbr=nbr)
    level += 1

    # _s_a, _s_b are just temporary place holders
    bottom = True
    scores = []
    for cid in cids:
        if cid in bfs:
            bottom = False
            _score = score(read, bfs[cid])
            if _score > score_cutoff:
                scores.append(_score)
            else:
                scores.append(None)

    cid_scores = zip(cids, scores)
    if DEBUG:
        logging.debug('score_cutoff: {cutoff}, current_bf: {bf_id}, level: {level}, '
                      '{cid_scores}'.format(level=level, bf_id=bf_id,
                                            cutoff=score_cutoff,
                                            cid_scores=cid_scores))

    cid_scores = [__ for __ in cid_scores if __[1] is not None]
    # Don't change the order of if in the following code! They must be this
    # way, or use elif with less eligibility. And don't forget to return after
    # each if

    # bf_id is the bottom of the bfs tree or score is too low for all children bfs
    if not cid_scores:
        hit.append(bf_id)
        if not bottom:
            with open(NO_MATCH, 'ab') as opf:
                opf.write('{0}\n'.format(read))
        return

    # find cids with the SAME MAX scores
    max_score = 0
    max_cids = []               # cids with max_scores
    for (c, s) in cid_scores:
        if s >= max_score:
            if s > max_score:
                max_score = s
                max_cids = [c]
            else:
                max_cids.append(c)

    for __ in max_cids:
        get_bf(read, bfs, nbr, score_cutoff, hit, __, level)
        return


def worker(pid, queue, bfs, nbr, score_cutoff, res_count, lock):
    score_cutoff = score_cutoff.value
    res = {}                    # hold result for this process

    counter = 0
    while True:
        read = queue.get()
        hit = []

        get_bf(read, bfs, nbr, score_cutoff, hit)

        for bf_id in hit:
            if bf_id in res:
                res[bf_id] += 1
            else:
                res[bf_id] = 1

        counter += 1
        if DEBUG:
            interval = 100
        else:
            interval = 5000
        if counter % interval == 0:
            logging.info('pid: {0}, analyzed {1} reads'.format(pid, counter))

        if queue.empty():
            more_task = False   # indicate whether there are more tasks or not
            if DEBUG:
                n_loop = 1
            else:
                n_loop = 10
            for i in xrange(1, n_loop + 1):
                logging.info('pid: {0}: waiting {1}s for queue to be non-empty'.format(pid, i))
                time.sleep(i)
                if not queue.empty():
                    logging.info('pid: {0}: more tasks found after {1}s '
                                 'wait, back to work'.format(pid, (i + 1) * i / 2))
                    more_task = True
                    break

            if more_task:
                continue
            else:
                total_seconds = (n_loop + 1) * n_loop / 2.
                logging.info('pid: {0}: still empty after {1}s wait for '
                             'the queue, break the while loop'.format(
                                 pid, total_seconds))
                break

    with lock:
        for k in res.keys():
            if k in res_count.keys():
                res_count[k] += res[k]
            else:
                res_count[k] = res[k]


def fetch_reads(*fq_gzs):
    for fq_gz in fq_gzs:
        logging.info('working on {0}'.format(fq_gz))
        with gzip.open(fq_gz) as inf:
            for rec in SeqIO.parse(inf, "fastq"):
                yield str(rec.seq)


if __name__ == "__main__":
    if os.path.exists(NO_MATCH):
        os.remove(NO_MATCH)

    if DEBUG:
        db_file = 'debug_db/combined.db'
    else:
        db_file = 'db/combined.db'

    # init shared variables among multiple processes
    # https://docs.python.org/2/library/multiprocessing.html#sharing-state-between-processes
    score_cutoff = Value('d', SCORE_CUTOFF) # value is faster than manager
    manager = Manager()
    # bfs = manager.dict(load_bfs(db_file))
    bfs = load_bfs(db_file)
    res_count = manager.dict() # to hold results of a dictionary of (bf_id: count)
    lock = Lock()

    # init queue for input data 
    queue = Queue()

    # start multiple processes
    procs = []
    for pid in range(NUM_CPUS - 1):
        proc = Process(target=worker, args=(pid, queue, bfs, NBR, score_cutoff,
                                            res_count, lock))
        proc.daemon = True
        procs.append(proc)
        proc.start()

    if DEBUG:
        from test_data import TEST_READS
        input_reads = TEST_READS
    else:
        fq_gzs = [
            '/projects/btl2/zxue/microorganism_profiling/real/BALC7.Clinical/reads_1.fq.gz',
            '/projects/btl2/zxue/microorganism_profiling/real/BALC7.Clinical/reads_2.fq.gz'
        ]
        input_reads = fetch_reads(*fq_gzs)

    # enqueue input data
    freq, next_freq = 1, 10
    for k, read in enumerate(input_reads):
        k_plus_1 = k + 1
        if k_plus_1 <= next_freq:
            if k_plus_1 % freq == 0:
                logging.info('enqueuing {0}th read'.format(k_plus_1))
                logging.info(read)
        else:
            freq *= 10
            next_freq = freq * 10            
        queue.put(read)
    logging.info('enqueued {0} reads in total'.format(k_plus_1))

    for proc in procs:
        proc.join()

    df = pd.DataFrame.from_records(res_count.items(), columns=['bf_id', 'hit_count'])
    df.sort('hit_count', ascending=False, inplace=True)
    df.to_csv('res_count.csv', index=False)

