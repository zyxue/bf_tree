from __future__ import unicode_literals

import bitarray
import sqlite3


from objs import BloomFilter
import utils as U

import logging
logging.basicConfig(level=logging.DEBUG,
                    # filename='app.log',
                    # filemode='w',
                    format='%(asctime)s|%(levelname)s|%(message)s')

K_MER_SIZE = 20

# seq_id: 0
# TEST_READ = 'ATTCTTTCGCTTAACTCTGTGTTTGCGTTAAACCTTCTGTGGTAGGTGTCTAGTCTACAGGCCTAACTGTTACTTGATTTCTTGTCAAGTCAGGCAGTGGAGTGACGACAAGCGATTGTTGCCGCTAGACTTAAAGATTATTCACCTCGA'
# mutated seq_id: 0
# TEST_READ = 'ATTCTTTCGCTTAACTCTGTGTTTGCGTTAAACCTTCTGTGGTAGGGGTCTAGTCTACAGGCCTAACTGTTACTTGATTTCTTGTCAAGTCAGGCAGTGGAGTGACGACAAGCGATTGTAGCCGCTAGACTTAAAGATTATTCACCTCGA'
# seq_id: 8
# TEST_READ = 'AGAACGCAGCGAAGCGCGAAATGTAGTGTGAATCGCAGAACATTGTGAATCATCGAATCTTTGAACGCACATTGCGCCTCCCGTTACTGGGAGGCATGCCTGTCTGAGCGTCATTTAAAACGTAGCACCCCACCAAGGGTCTGTCTTGGG'
# seq_id: 88
TEST_READ = 'CTCGCATCGATGAAGAACGCAGCGAAGCGCGAAATGTAATGTGAATCGCAGAACATTGTGAATCATCGAATCTTTGAACGCACATTGCGCCTCCAGTAACTCGGAGGCATGCCTGTCTGAGCGTCATTTAAAACCATGGACCCGACCGGG'


SCORE_CUTOFF = 0.5

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
    
    bfs = {}
    for bf_id, size, hash_count, bf in cur.fetchall():
        bfs[bf_id] = BloomFilter(size, hash_count, bf)
    conn.close()

    logging.info('loading is done'.format(db_file))
    return bfs


COMPLEMENT_DD = {
    'A': 'T',
    'C': 'G',
    'T': 'A',
    'G': 'C',
    'N': 'N'
}
def reverse_complement(seq):
    return ''.join([COMPLEMENT_DD[x] for x in seq[::-1]])



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
    f_read = TEST_READ
    r_read = reverse_complement(f_read)

    f_s = calc_score(f_read, K_MER_SIZE, bf)
    r_s = calc_score(r_read, K_MER_SIZE, bf)

    # print f_s, r_s
    return max(f_s, r_s)


def main():
    db_file = 'lele_db/combined.db'
    # db_file = 'db.bk/2.db'

    bfs = load_bfs(db_file)

    # from pprint import pprint
    # pprint(bfs)

    res = get_bf(TEST_READ, bfs, SCORE_CUTOFF)
    print 'bf_id: ', res


    # bf = bitarray.bitarray()
    # bf.frombytes(str(bfs[1]))
    # print bf


def get_bf(read, bfs, score_cutoff, bf_id=0, level=0):
    """if calculated score < score_cutoff, then stop searching"""
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

    logging.debug('{a}: {s_a}, {b}: {s_b}, level: {level}, '
                  'score_cutoff: {cutoff}'.format(a=a, b=b, s_a=_s_a, s_b=_s_b,
                                                  level=level, cutoff=score_cutoff))

    # Don't change the order of if in the following code! They must be this
    # way, or use elif with less eligibility

    # bf_id is the bottom of the bfs tree or score is two low for both children bfs
    if s_a is None and s_b is None:
        return bf_id

    if s_a is not None and s_b is not None:
        if s_a > s_b:
            return get_bf(read, bfs, score_cutoff, a, level)
        elif s_a < s_b:
            return get_bf(read, bfs, score_cutoff, b, level)
        else:
            # multimatch: pass through both channels
            get_bf(read, bfs, score_cutoff, a, level)
            get_bf(read, bfs, score_cutoff, b, level)
            return 

    if s_a is not None:
        if s_a > score_cutoff:
            return get_bf(read, bfs, score_cutoff, a, level)
        else:
            return bf_id

    if s_b is not None:
        if s_b > score_cutoff:
            return get_bf(read, bfs, score_cutoff, b, level)
        else:
            return bf_id


main()
