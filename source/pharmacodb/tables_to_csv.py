from __future__ import print_function
import MySQLdb
import csv

user = 'root'  # your username
passwd = 'pass'  # your password
host = 'localhost'  # your host or localhost if running locally
db = 'pharmacodb_development'  # database where your table is stored
dest_folder = '/var/lib/mysql-files/'  # destination folder for the files to be written

con = MySQLdb.connect(user=user, passwd=passwd, host=host, db=db)
cursor = con.cursor()

tables = ["cell_tissues", "cellosaurus", "cells", "dataset_cells", "datasets", "dose_responses", "drug_annots",
          "drugs", "experiments", "profiles", "source_cell_names", "source_drug_names", "source_statistics",
          "source_tissue_names", "sources", "tissues"]
for table in tables:
    print("dumping table:", table)
    with open(dest_folder + table + '.csv', 'w') as f:
        writer = csv.writer(f)

        query = "SELECT COLUMN_NAME FROM information_schema.COLUMNS C WHERE table_name = '%s';" % table
        cursor.execute(query)
        columns = []
        for c in cursor.fetchall():
            columns.append(c[0])
        writer.writerow(columns)

        query = "SELECT * FROM %s;" % table
        cursor.execute(query)
        for row in cursor.fetchall():
            writer.writerow(row)
