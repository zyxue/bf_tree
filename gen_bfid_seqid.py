import pandas as pd
import sqlite3

conn = sqlite3.connect('db/v4/nbr8/combined.db')
cursor = conn.cursor()
cursor.execute('select bf_id, seq_id from bloomfilter where seq_id IS NOT NULL;')

res = cursor.fetchall() 

bf_df = pd.DataFrame.from_records(res, columns=['bf_id', 'seq_id'])
seq_df = pd.read_csv('/projects/btl2/zxue/microorganism_profiling/libraries_assesment/comparison/v4/silva_gg_rdp_unite_no_eukaryotes.csv.gz', index_col=0, compression='gzip')
seq_df.drop(['seq'], axis=1, inplace=True)
seq_df = seq_df.merge(bf_df, on='seq_id')

# seq_df = seq_df.merge(bf_seq_df, on='seq_id')[['bf_id', 'seq_id', 'seq_len', 'order', 'family', 'genus', 'species']]
