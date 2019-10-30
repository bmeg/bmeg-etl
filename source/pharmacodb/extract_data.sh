#!/bin/bash
set -e

mysqld &

echo "waiting for mysql to start..."
while ! mysqladmin ping --silent; do
    sleep 1
done
echo "done."

echo "extracting tables from pharmacodb_development..."
apt-get install -y python-mysqldb
python tables_to_csv.py
echo "done."

chmod 666 *.csv
