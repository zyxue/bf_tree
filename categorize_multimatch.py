from __future__ import unicode_literals

import os
import csv
import time
import sqlite3
import random
import pandas as pd
from multiprocessing import Manager, Value, Queue, Process, Lock
import gzip
import logging

from Bio import SeqIO

from objs import BloomFilter
import utils as U
from constants import COMPLEMENT_DD

from settings import DEBUG

if DEBUG:
    logging.basicConfig(level=logging.DEBUG,
                        format='%(asctime)s|%(levelname)s|%(message)s')
else:
    logging.basicConfig(level=logging.DEBUG,
                        filename='categorize.log', filemode='w',
                        format='%(asctime)s|%(levelname)s|%(message)s')


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


def score(read, k_mer_size, bf):
    """
    read: should be a string, not SeqRecord
    """
    f_read = read
    r_read = reverse_complement(f_read)

    f_s = calc_score(f_read, k_mer_size, bf)
    r_s = calc_score(r_read, k_mer_size, bf)

    # print f_s, r_s
    return max(f_s, r_s)


def get_bf(read, bfs, nbr, k_mer_size, score_cutoff, hit, no_match_output,
           bf_id=0, bf_score=0, level=0):
    """
    if calculated score < score_cutoff, then stop searching

    :param hit: a list of tuples like (bf_id, score) storing result for this read
    """

    cids = U.calc_children_id(bf_id, level, nbr=nbr)
    level += 1

    # _s_a, _s_b are just temporary place holders
    bottom = True
    scores = []
    for cid in cids:
        if cid in bfs:
            bottom = False
            _score = score(str(read.seq), k_mer_size, bfs[cid])
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
        if bottom:
            hit.append((bf_id, bf_score))
        # else:
        #     with gzip.open(no_match_output, 'ab') as opf:
        #         opf.write(read.format('fastq'))
        return

    for (bf_id, bf_score) in cid_scores:
        get_bf(read, bfs, nbr, k_mer_size, score_cutoff, hit, no_match_output,
               bf_id, bf_score, level)
    return


def output_res(res, output_file):
    with gzip.open(output_file, 'ab') as opf:
        csvwriter = csv.writer(opf)
        csvwriter.writerows(res)


def worker(pid, queue, lock, bfs, nbr, k_mer_size, score_cutoff, prefix, output_dir):
    score_cutoff = score_cutoff.value
    score_output = os.path.join(output_dir, '{pid}_{prefix}_scores.csv.gz'.format(**locals()))
    no_match_output = os.path.join(output_dir, '{pid}_{prefix}_no_match.fq.gz'.format(**locals()))

    for _ in [score_output, no_match_output]:
        if os.path.exists(_):
            os.remove(_)
                
    res = []      # hold list of (read_id, bf_id, score) as itermediate results

    counter = 0
    while True:
        read = queue.get()
        hit = []

        get_bf(read, bfs, nbr, k_mer_size, score_cutoff, hit, no_match_output)

        for (bf_id, score) in hit:
            res.append((read.read_id, bf_id, score))

        if DEBUG:
            write_trigger = 1
        else:
            write_trigger = random.randint(10000, 90000)
        # to prevent multiple processes from writing to the disk simultaneously
        if len(res) > write_trigger:
            output_res(res, score_output)
            res = []

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

                output_res(res, score_output)
                break



def fetch_reads(*fq_gzs):
    count = 0
    for fq_gz in fq_gzs:
        logging.info('working on {0}'.format(fq_gz))
        with gzip.open(fq_gz) as inf:
            for rec in SeqIO.parse(inf, "fastq"):
                rec.read_id = count
                count += 1
                # yield str(rec.seq)
                yield rec


if __name__ == "__main__":
    # inputs include: 1, 2, 3 as shown belown

    # 1. input fq files, output_dir (default to the same dir of fq files),
    # prefix (drived from fq files or arbitrarily set), output file names
    # (defined based on prefix and output_dir)
    if DEBUG:
        fq_gzs = ['test_data.fq.gz']
        output_dir = 'debug_results/v4'
    else:
        fq_gzs = [
            '/projects/btl2/zxue/microorganism_profiling/real/BALC7.Clinical/reads_1.fq.gz',
            '/projects/btl2/zxue/microorganism_profiling/real/BALC7.Clinical/reads_2.fq.gz'

            # '/projects/btl2/zxue/microorganism_profiling/real/minnow/SRR1561864_1.fq.gz',
            # '/projects/btl2/zxue/microorganism_profiling/real/minnow/SRR1561864_2.fq.gz'
        ]
        output_dir = 'results/v4'
    prefix = 'BALC7'
    input_reads = fetch_reads(*fq_gzs)


    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    output_no_match = os.path.join(output_dir, '{0}_no_match.fq'.format(prefix))


    # 2. database of bloomfilters
    if DEBUG:
        db_file = 'debug_db/v4/nbr8/combined.db'
    else:
        db_file = 'db/v4/nbr8/combined.db'

    # 3. configuration parameters
    if DEBUG:
        num_cpus = 4
    else:
        num_cpus = 32               # given 32 cores
    
    score_cutoff = 0.5          # future: get from config
    k_mer_size = 20             # future: get from db
    nbr = 8                     # future: get from db
    
    
    # init shared variables among multiple processes
    # https://docs.python.org/2/library/multiprocessing.html#sharing-state-between-processes
    score_cutoff = Value('d', score_cutoff) # value is faster than manager
    manager = Manager()
    # bfs = manager.dict(load_bfs(db_file))
    bfs = load_bfs(db_file)
    lock = Lock()

    # init queue for input data 
    queue = Queue()

    # start multiple processes
    procs = []
    for pid in range(num_cpus - 1):
        proc = Process(target=worker, args=(pid, queue, lock, bfs, nbr, k_mer_size,
                                            score_cutoff, prefix, output_dir))
        proc.daemon = True
        procs.append(proc)
        proc.start()

    # enqueue input data
    freq, next_freq = 1, 10
    for k, read in enumerate(input_reads):
        k_plus_1 = k + 1
        if k_plus_1 <= next_freq:
            if k_plus_1 % freq == 0:
                logging.info('enqueuing {0}th read'.format(k_plus_1))
        else:
            freq *= 10
            next_freq = freq * 10            
        queue.put(read)
    logging.info('enqueued {0} reads in total'.format(k_plus_1))

    for proc in procs:
        proc.join()

    # df = pd.DataFrame.from_records(res_count.items(), columns=['bf_id', 'hit_count'])
    # df.sort('hit_count', ascending=False, inplace=True)
    # df.to_csv('results/v4/BALC7.csv', index=False)
    # df.to_csv('results/v4/SRR1561864.csv', index=False)
