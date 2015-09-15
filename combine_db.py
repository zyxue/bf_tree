import os
import glob
import sqlite3

db_dir = 'lele_db'

output_db = os.path.join(db_dir, 'combined.db')

if os.path.exists(output_db):
    os.remove(output_db)


dbs = glob.glob(os.path.join(db_dir, '*'))

conn = sqlite3.connect(output_db)
cursor = conn.cursor()

# cursor.execute("CREATE TABLE bloomfilter (bf_id INTEGER PRIMARY KEY, num_uniq_kmers INTEGER, size INTEGER, hash_count INTEGER, fpr REAL, bitarray BLOB, seq_id)")
cursor.execute("CREATE TABLE bloomfilter (bf_id INTEGER PRIMARY KEY, bitarray BLOB, seq_id)")


for db in dbs:
    print('working on {0}'.format(db))
    cursor.execute('attach "{0}" as toCombine'.format(db))
    cursor.execute('insert into bloomfilter select * from toCombine.bloomfilter'.format(db))
    cursor.execute('detach toCombine'.format(db))

conn.commit()
