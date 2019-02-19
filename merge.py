import sqlite3
import pandas as pd


def gen_bf_df(bf_db):
    conn = sqlite3.connect(bf_db)
    cursor = conn.cursor()
    cursor.execute('select bf_id, seq_id from bloomfilter where seq_id IS NOT NULL;')

    res = cursor.fetchall() 
    
    bf_df = pd.DataFrame.from_records(res, columns=['bf_id', 'seq_id'])
    return bf_df
    

# hit_count = 'results/v1/BALC7.csv'
# bf_db = 'db/v1/nbr8/combined.db'
# seq_db = '/projects/btl2/zxue/microorganism_profiling/libraries_assesment/gg_unite/gg_unite.csv.gz'

# df = pd.read_csv(hit_count)
# bf_df = gen_bf_df(bf_db)

# seq_df = pd.read_csv(seq_db, compression='gzip', index_col=0)
# seq_df['seq_id'] = seq_df.index.values

# merged = df.merge(bf_df, on='bf_id').merge(seq_df, on='seq_id')



hit_count = 'results/v2/BALC7.csv'
bf_db = 'db/v2/nbr8/combined.db'
seq_db = '/projects/btl2/zxue/microorganism_profiling/libraries_assesment/comparison/genbank_db/filtered_gb_short_seq.csv.gz'
lineage_db = '/projects/btl2/zxue/microorganism_profiling/libraries_assesment/ncbi_taxonomy/species_lineages/ncbi_lineages.csv.gz'

hit_df = pd.read_csv(hit_count)
bf_df = gen_bf_df(bf_db)
seq_df = pd.read_csv(seq_db, compression='gzip')
seq_df['seq_id'] = seq_df.index.values
lineage_df = pd.read_csv(lineage_db)

merged = hit_df.merge(bf_df, on='bf_id').merge(seq_df, on='seq_id').merge(lineage_df, on='tax_id')
