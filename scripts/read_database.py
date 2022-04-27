from dataclasses import replace
import sqlite3
import pandas as pd


con = sqlite3.connect("db_test.db")
#df = pd.read_csv("2022-04-07-desh-meta.tsv", sep="\t")
#df.to_sql("metadata", con=con, if_exists="replace", index=False)
cur = con.cursor()
#cur.execute("ALTER TABLE metadata ADD dNdS real")
#con.commit()
for row in cur.execute("SELECT dNdS from metadata WHERE dNdS <> 'None' ORDER BY SAMPLING_DATE"):
    print(row)
con.close()