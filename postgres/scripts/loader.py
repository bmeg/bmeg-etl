"""
loads json vertex files into postgres table edge(gid, label, from, to, data)
"""

import dataset
import ujson
import logging
import yaml
import Queue
import threading

# DDL
DROP_TABLES = """
DROP TABLE IF EXISTS vertex;
DROP TABLE IF EXISTS edge;
"""

CREATE_TABLES = """
CREATE TABLE IF NOT EXISTS  vertex (
 gid varchar not null,
 label varchar not null,
 data jsonb
);

CREATE TABLE IF NOT EXISTS  edge (
 gid varchar not null,
 label varchar not null,
 "from" varchar not null,
 "to" varchar not null,
 data jsonb
);
"""

CREATE_INDEXES = """
CREATE INDEX vertex_gid ON vertex (gid);
CREATE INDEX edge_label_from_to ON edge (label, "from", "to");
CREATE INDEX edge_label_to_from ON edge (label, "from", "to");
ANALYZE vertex ;
ANALYZE edge ;
"""


# list of files for import
EDGE_FILES = """
outputs/gtex/gtex.AliquotFor.Edge.json
outputs/gtex/gtex.BiosampleFor.Edge.json
outputs/gtex/gtex.ExpressionOf.Edge.json
outputs/gtex/gtex.InProject.Edge.json
""".strip().split()

VERTEX_FILES = """
outputs/gtex/gtex.Aliquot.Vertex.json
outputs/gtex/gtex.Biosample.Vertex.json
outputs/gtex/gtex.Expression.Vertex.json
outputs/gtex/gtex.Individual.Vertex.json
outputs/gtex/gtex.Project.Vertex.json
""".strip().split()

# log setup
logging.getLogger().setLevel(logging.INFO)


# connection for dataset high level connection
# https://dataset.readthedocs.io/en/latest/
def construct_pg_url(user, password, host, port, database):
    if not password:
        return "postgresql://" + user + '@' + host + ':' + str(port) + '/' + database
    return "postgresql://" + user + ":" + password + '@' + host + ':' + str(port) + '/' + database


with open("postgres/scripts/config.yml", 'r') as stream:
    config = yaml.load(stream)

pgconn = dataset.Database(url=construct_pg_url())


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
        with open(f) as ins:
            for line in ins:
                c += 1
                t += 1
                obj = ujson.loads(line)
                for k in keys_to_delete:
                    del obj[k]
                yield transform(f, obj)
                if c % batch_size == 0:
                    c = 0
                    logging.info('loaded {}'.format(t))


logging.info('(re) creating tables')
execute(pgconn, [config.ddl])

# connect to table, load it and log count

SENTINEL = 'XXXXXX'


def writer_worker(q, table_name):
    """ write to table name """
    t = pgconn[table_name]
    t.insert_many(iter(q.get, SENTINEL), chunk_size=100)


def reader_worker(q, files):
    """ write to q from files """
    for row in rows(files):
        q.put(row)
    q.put(SENTINEL)


# create queues and threads
vertex_q = Queue(maxsize=100)
edge_q = Queue(maxsize=100)
threads = []

for i in range(10):
    t = threading.Thread(target=writer_worker, args=(vertex_q, 'vertex'))
    t.start()
    threads.append(t)

for i in range(10):
    t = threading.Thread(target=writer_worker, args=(edge_q, 'edge'))
    t.start()
    threads.append(t)

t = threading.Thread(target=reader_worker, args=(vertex_q, config.vertex_files))
t.start()
threads.append(t)
t = threading.Thread(target=reader_worker, args=(edge_q, config.edge_files))
t.start()
threads.append(t)

# Blocks until all items in the queue have been gotten and processed.
for q in [vertex_q, edge_q]:
    q.join()
# block until all tasks are done
for t in threads:
    t.join()


logging.info('creating indexes')
execute(pgconn, [config.indexes])
logging.info('indexes created')
