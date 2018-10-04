"""
loads json vertex files into postgres table edge(gid, label, from, to, data)
"""

import dataset
import ujson
import logging
import yaml
from queue import Queue
import threading
import os
import types
from bmeg.ioutils import reader

# log setup
logging.getLogger().setLevel(logging.INFO)


# connection for dataset high level connection
# https://dataset.readthedocs.io/en/latest/
def construct_pg_url(user, host, port, database, password=None):
    if not password:
        return "postgresql://" + user + '@' + host + ':' + str(port) + '/' + database
    return "postgresql://" + user + ":" + password + '@' + host + ':' + str(port) + '/' + database


with open("{}/config.yml".format(os.path.dirname(os.path.realpath(__file__))), 'r') as stream:
    config = yaml.load(stream)

config = types.SimpleNamespace(**config)
config.edge_files = config.edge_files.strip().split()
config.vertex_files = config.vertex_files.strip().split()

pgconn = dataset.Database(url=construct_pg_url(**config.postgres), engine_kwargs={'pool_size': 30, 'max_overflow': 20})


def execute(pgconn, commands):
    """ exec statement on connection and commit """
    pgconn.begin()
    for command in commands:
        pgconn.query(command)
    pgconn.commit()


def gtex_expression_data(row):
    """
    flatten gtex expression data from data.values.ENSG00000223972 = 0.1081
    to data.values = [ {name:'ENSG00000223972', value:0.1081}, ...]
    """
    # row['data']['values'] = [[k,row['data']['values'][k]] for k in row['data']['values'].keys()]
    row['data']['values'] = [{'name': k, 'value': row['data']['values'][k]} for k in row['data']['values'].keys()]
    return row


def transform(file, row):
    """ transform if necessary """
    # if file == 'outputs/gtex/gtex.Expression.Vertex.json':
    #     return gtex_expression_data(row)
    return row


def rows(files, keys_to_delete=['_id'], batch_size=1000):
    """
    generator: read in all rows from all files,
    remove keys_to_delete before yielding
    log message every batch_size
    """
    t = c = 0
    for f in files:
        logging.info('reading {}'.format(f))
        t = 0
        try:
            with reader(f) as ins:
                for line in ins:
                    c += 1
                    t += 1
                    obj = ujson.loads(line)
                    for k in keys_to_delete:
                        del obj[k]
                    yield transform(f, obj)
                    if c % batch_size == 0:
                        c = 0
                        logging.info('loaded {} {}'.format(f, t))
        except Exception as e:
            logging.exception(e)
            logging.error(f)


logging.info('(re) creating tables')
execute(pgconn, [config.ddl])

# connect to table, load it and log count

SENTINEL = 'XXXXXX'


def writer_worker(q, table_name):
    """ write to table name """
    logging.info('writer worker started')
    t = pgconn[table_name]
    t.insert_many(iter(q.get, SENTINEL), chunk_size=1000)
    logging.info('writer worker done')


def reader_worker(q, files):
    """ write to q from files """
    for row in rows(files):
        q.put(row)


# create queues and threads
vertex_q = Queue(maxsize=2000000)
edge_q = Queue(maxsize=2000000)
writer_threads = []
reader_threads = []

for i in range(10):
    t = threading.Thread(target=writer_worker, args=(vertex_q, 'vertex'))
    t.start()
    writer_threads.append(t)

for i in range(10):
    t = threading.Thread(target=writer_worker, args=(edge_q, 'edge'))
    t.start()
    writer_threads.append(t)

for vertex_file in config.vertex_files:
    t = threading.Thread(target=reader_worker, args=(vertex_q, [vertex_file]))
    t.start()
    reader_threads.append(t)

for edge_file in config.edge_files:
    t = threading.Thread(target=reader_worker, args=(edge_q, [edge_file]))
    t.start()
    reader_threads.append(t)

# Blocks until all items in the queue have been gotten and processed.
logging.info('waiting on queues')
for q in [vertex_q, edge_q]:
    q.join()
# block until all tasks are done
logging.info('waiting on reader_threads')
for t in reader_threads:
    t.join()
# block until all reader tasks are done
for t in reader_threads:
    # tell writers to exit
    for q in [vertex_q, edge_q]:
        q.put(SENTINEL)
# block until all writer tasks are done
logging.info('waiting on writer_threads')
for t in writer_threads:
    t.join()


logging.info('creating indexes')
execute(pgconn, [config.indexes])
logging.info('indexes created')
