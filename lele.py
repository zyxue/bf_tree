from bitarray import bitarray
import sqlite3

a  = bitarray('1001011')

conn = sqlite3.connect(":memory:")
# conn = sqlite3.connect('tmp.db')
cursor = conn.cursor()

cursor.execute("CREATE TABLE t (id INTEGER, x BLOB)")

cursor.execute("INSERT INTO t values (?, ?)", (1, sqlite3.Binary(a.tobytes())))
cursor.execute("select * from t")

rec = cursor.fetchone()

conn.commit()
conn.close()

b = bitarray.frombytes(rec[1])
