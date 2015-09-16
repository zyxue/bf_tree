from __future__ import unicode_literals

import os
import time
import bitarray
import sqlite3
import pandas as pd
from multiprocessing import Manager, Value, Queue, Process, Lock


from objs import BloomFilter
import utils as U
from constants import COMPLEMENT_DD

import logging
logging.basicConfig(level=logging.DEBUG,
                    # filename='app.log',
                    # filemode='w',
                    format='%(asctime)s|%(levelname)s|%(message)s')

from settings import DEBUG


K_MER_SIZE = 20
SCORE_CUTOFF = 0.5

if DEBUG:
    NUM_CPUS = 4
    DB_OUTPUT_DIR = 'debug_db'
else:
    NUM_CPUS = 32         # given 32 cores
    DB_OUTPUT_DIR = 'db'
if not os.path.exists(DB_OUTPUT_DIR):
    os.mkdir(DB_OUTPUT_DIR)


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
        k_plus_1 = k + 1
        if k_plus_1 < next_freq:
            if k_plus_1 % freq == 0:
                logging.info('working on {0}th bf'.format(k_plus_1))
        else:
            freq *= 10
            next_freq = freq * 10            
        bfs[bf_id] = BloomFilter(size, hash_count, bf)
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


def get_bf(read, bfs, score_cutoff, hit, bf_id=0, level=0):
    """
    if calculated score < score_cutoff, then stop searching

    :param hit: a list that whole hit bf_ids
    """

    a, b = U.calc_children_id(bf_id, level)
    level += 1

    # _s_a, _s_b are just temporary place holders
    _s_a, _s_b, s_a, s_b = None, None, None, None
    if a in bfs:
        _s_a = score(read, bfs[a])
        if _s_a >= score_cutoff:
            s_a = _s_a
    if b in bfs:
        _s_b = score(read, bfs[b])
        if _s_b >= score_cutoff:
            s_b = _s_b

    # logging.debug('score_cutoff: {cutoff}, level: {level}, '
    #               '({a}): {s_a}, ({b}): {s_b}'.format(
    #                   a=a, b=b, s_a=_s_a, s_b=_s_b,
    #                   level=level, cutoff=score_cutoff))

    # Don't change the order of if in the following code! They must be this
    # way, or use elif with less eligibility. And don't forget to return after
    # each if

    # bf_id is the bottom of the bfs tree or score is two low for both children bfs
    if s_a is None and s_b is None:
        hit.append(bf_id)
        return

    if s_a is not None and s_b is not None:
        if s_a > s_b:
            get_bf(read, bfs, score_cutoff, hit, a, level)
            return
        elif s_a < s_b:
            get_bf(read, bfs, score_cutoff, hit, b, level)
            return
        else:
            # multimatch: pass through both channels
            get_bf(read, bfs, score_cutoff, hit, a, level)
            get_bf(read, bfs, score_cutoff, hit, b, level)
            return

    # Compared to the above block, this block encourages multiMatch, which may
    # increase the complexity of results interpretation

    # if s_a is not None and s_b is not None:
    #     # multimatch: pass through both channels
    #     get_bf(read, bfs, score_cutoff, hit, a, level)
    #     get_bf(read, bfs, score_cutoff, hit, b, level)
    #     return 

    if s_a is not None:
        if s_a > score_cutoff:
            get_bf(read, bfs, score_cutoff, hit, a, level)
        else:
            hit.append(bf_id)
        return

    if s_b is not None:
        if s_b > score_cutoff:
            get_bf(read, bfs, score_cutoff, hit, b, level)
        else:
            hit.append(bf_id)
        return


def worker(pid, queue, bfs, score_cutoff, res_count, lock):
    score_cutoff = score_cutoff.value
    res = {}                    # hold result for this process
    while True:
        read = queue.get()
        hit = []
        get_bf(read, bfs, score_cutoff, hit)
        for bf_id in hit:
            if bf_id in res:
                res[bf_id] += 1
            else:
                res[bf_id] = 1

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


if __name__ == "__main__":
    if DEBUG:
        db_file = 'debug_db/combined.db'
    else:
        db_file = 'db/combined.db'


    # init shared variables among multiple processes
    # https://docs.python.org/2/library/multiprocessing.html#sharing-state-between-processes
    score_cutoff = Value('d', SCORE_CUTOFF) # value is faster than manager
    manager = Manager()
    bfs = manager.dict(load_bfs(db_file))
    res_count = manager.dict() # to hold results of a dictionary of (bf_id: count)
    lock = Lock()

    # init queue for input data 
    queue = Queue()

    # start multiple processes
    procs = []
    for pid in range(NUM_CPUS - 1):
        proc = Process(target=worker, args=(pid, queue, bfs, score_cutoff,
                                            res_count, lock))
        proc.daemon = True
        procs.append(proc)
        proc.start()

    if DEBUG:
        from pprint import pprint
        pprint(bfs)

    from test_data import TEST_READS, TEST_READ_4

    # enqueue input data
    freq, next_freq = 1, 10
    for k, read in enumerate(TEST_READS):
        k_plus_1 = k + 1
        if k_plus_1 < next_freq:
            if k_plus_1 % freq == 0:
                logging.info('enqueuing {0}th read'.format(k_plus_1))
                logging.info(read)
        else:
            freq *= 10
            next_freq = freq * 10            
        queue.put(read)

    for proc in procs:
        proc.join()
            
    for k in sorted(res_count.items(), key=lambda x: x[1], reverse=True):
        print k

    df = pd.DataFrame.from_records(res_count.items(), columns=['bf_id', 'seq_id'])
    df.sort('seq_id', ascending=False, inplace=True)
    df.to_csv('res_count.csv', index=False)

