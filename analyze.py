from Bio import SeqIO
import gzip
import logging
logging.basicConfig(level=logging.INFO,
                    # filename='app.log',
                    # filemode='w',
                    format='%(asctime)s|%(levelname)s|%(message)s')


res = {}

infile1 = '/projects/btl2/zxue/microorganism_profiling/real/BALC7.Clinical/reads_1.fq.gz'
logging.info('working on {0}'.format(infile1))
with gzip.open(infile1) as inf:
    for k, rec in enumerate(SeqIO.parse(inf, "fastq")):
        if (k + 1) % 100 == 0:
            items = res.items()
            for ii in sorted(items, key=lambda x: x[1], reverse=True)[:20]:
                print ii
            logging.info('working on {0}th read'.format(k + 1))
        # print __
        hit = []
        for (bf_id, level) in [(3, 2), (4, 2), (5, 2), (6, 2)]:
            get_bf(str(rec.seq), bfs, 0.5, hit, bf_id, level)
        for bf_id in hit:
            if bf_id in res:
                res[bf_id] += 1
            else:
                res[bf_id] = 1

    
infile2 = '/projects/btl2/zxue/microorganism_profiling/real/BALC7.Clinical/reads_2.fq.gz'
logging.info('working on {0}'.format(infile2))
with gzip.open(infile2) as inf:
    for k, rec in enumerate(SeqIO.parse(inf, "fastq")):
        if (k + 1) % 10000 == 0:
            logging.info('working on {0}th read'.format(k + 1))
        # print __
        hit = []
        for (bf_id, level) in [(3, 2), (4, 2), (5, 2), (6, 2)]:
            get_bf(str(rec.seq), bfs, 0.5, hit, bf_id, level)
        for bf_id in hit:
            if bf_id in res:
                res[bf_id] += 1
            else:
                res[bf_id] = 1

