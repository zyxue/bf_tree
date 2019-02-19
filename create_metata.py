import sqlite3


output_db = '/projects/btl2/zxue/bf_tree/db/v2/nbr8/combined.db'

conn = sqlite3.connect(output_db)
cursor = conn.cursor()

# cursor.execute(SQL_CREATE_TABLE)

cursor.execute("INSERT INTO metadata values (?, ?, ?)", (K_MER_SIZE, NBR, FALSE_POSITIVE_RATE))

conn.commit()
