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
import argparse

# log setup
logging.getLogger().setLevel(logging.INFO)


# connection for dataset high level connection
# https://dataset.readthedocs.io/en/latest/
def construct_pg_url(user, host, port, database, password=None):
    if not password:
        return "postgresql://" + user + '@' + host + ':' + str(port) + '/' + database
    return "postgresql://" + user + ":" + password + '@' + host + ':' + str(port) + '/' + database


def execute(pgconn, commands):
    """ exec statement on connection and commit """
    pgconn.begin()
    for command in commands:
        pgconn.query(command)
    pgconn.commit()


def transform(file, row, keys_to_delete):
    """ transform if necessary """
    # if file == 'outputs/gtex/gtex.Expression.Vertex.json':
    #     return gtex_expression_data(row)
    for k in keys_to_delete:
        del row[k]
    return row


def rows(files, keys_to_delete=['_id'], batch_size=10000):
    """
    generator: read in all rows from all files,
    remove keys_to_delete before yielding
    log message every batch_size
    """
    t = c = 0
    for f in files:
        # logging.info('reading {}'.format(f))
        t = 0
        try:
            with reader(f) as ins:
                for line in ins:
                    c += 1
                    t += 1
                    obj = ujson.loads(line)
                    yield transform(f, obj, keys_to_delete)
                    if c % batch_size == 0:
                        c = 0
                        logging.info('loaded {} {}'.format(f, t))
            logging.info('loaded {}'.format(f))
        except Exception as e:
            logging.exception(e)
            logging.error(f)


def matrix_rows(files, keys_to_delete=['_id'], batch_size=10000):
    """
    generator: read in all rows from all files,
    remove keys_to_delete before yielding
    log message every batch_size
    """
    t = c = 0
    for f in files:
        # logging.info('reading {}'.format(f))
        t = 0
        try:
            with reader(f) as ins:
                for line in ins:
                    c += 1
                    t += 1
                    obj = ujson.loads(line)
                    for k in keys_to_delete:
                        del obj[k]
                    del obj['data']['values']
                    values = ujson.loads(line)['data']['values']
                    matrix = [{'gid': obj['gid'], 'name': k, 'value': values[k]} for k in values.keys() if values[k] != 0]
                    yield obj, matrix
                    if c % batch_size == 0:
                        c = 0
                        logging.info('loaded {} {}'.format(f, t))
            logging.info('loaded {}'.format(f))
        except Exception as e:
            logging.exception(e)
            logging.error(f)

# connect to table, load it and log count


SENTINEL = 'XXXXXX'


def writer_worker(pgconn, q, table_name):
    """ write to table name """
    logging.info('writer worker started')
    t = pgconn[table_name]
    t.insert_many(iter(q.get, SENTINEL), chunk_size=500)
    logging.info('writer worker done')


def reader_worker(q, files):
    """ write to q from files """
    for row in rows(files):
        q.put(row)


def matrix_reader_worker(vertex_q, matrix_q, files):
    """ write to q from files """
    for vertex, matrix in matrix_rows(files):
        vertex_q.put(vertex)
        matrix_q.put(vertex)


def main(dry, drop, index, config, workers=1):
    with open(config, 'r') as stream:
        config = yaml.load(stream)

    config = types.SimpleNamespace(**config)
    config.edge_files = config.edge_files.strip().split()
    config.vertex_files = config.vertex_files.strip().split()
    config.matrix_files = config.matrix_files.strip().split()
    all_files = config.vertex_files + config.edge_files + config.matrix_files

    pgconn = None
    if dry:
        logging.info('dry run: would have connected to {}'.format(construct_pg_url(**config.postgres)))
    else:
        pgconn = dataset.Database(url=construct_pg_url(**config.postgres), engine_kwargs={'pool_size': 30, 'max_overflow': 20})

    if drop:
        if dry:
            logging.info('dry run: would have dropped {}'.format([config.drop]))
        else:
            logging.info('dropping tables')
            execute(pgconn, [config.drop])
            logging.info('tables dropped')
    else:
        logging.info('skipping dropping tables')

    if dry:
        logging.info('dry run: would have created {}'.format([config.ddl]))
    else:
        logging.info('creating tables')
        execute(pgconn, [config.ddl])
        logging.info('tables created')

    for fname in all_files:
        assert os.path.isfile(fname), '{} does not exist'.format(fname)

    if dry:
        logging.info('dry run: read from {}'.format(all_files))
    else:
        # create queues and threads
        vertex_q = Queue(maxsize=2000)
        edge_q = Queue(maxsize=2000)
        matrix_q = Queue(maxsize=2000)
        all_queues = [vertex_q, edge_q, matrix_q]

        writer_threads = []
        reader_threads = []

        t = threading.Thread(target=writer_worker, args=(pgconn, vertex_q, 'vertex'))
        t.start()
        writer_threads.append(t)

        t = threading.Thread(target=writer_worker, args=(pgconn, edge_q, 'edge'))
        t.start()
        writer_threads.append(t)

        t = threading.Thread(target=writer_worker, args=(pgconn, matrix_q, 'matrix'))
        t.start()
        writer_threads.append(t)

        t = threading.Thread(target=reader_worker, args=(vertex_q, config.vertex_files))
        t.start()
        reader_threads.append(t)

        t = threading.Thread(target=reader_worker, args=(edge_q, config.edge_files))
        t.start()
        reader_threads.append(t)

        t = threading.Thread(target=matrix_reader_worker, args=(vertex_q, matrix_q, config.matrix_files))
        t.start()
        reader_threads.append(t)

        # Blocks until all items in the queue have been gotten and processed.
        logging.info('waiting on queues')
        for q in all_queues:
            q.join()
        # block until all tasks are done
        logging.info('waiting on reader_threads')
        for t in reader_threads:
            t.join()
        # block until all reader tasks are done
        for t in reader_threads:
            # tell writers to exit
            for q in all_queues:
                q.put(SENTINEL)
        # block until all writer tasks are done
        logging.info('waiting on writer_threads')
        for t in writer_threads:
            t.join()

    if index:
        if dry:
            logging.info('dry run: would have indexed {}'.format([config.indexes]))
        else:
            logging.info('creating indexes')
            execute(pgconn, [config.indexes])
            logging.info('indexes created')
    else:
        logging.info('skipping index')

    logging.info('Done!')


if __name__ == '__main__':  # pragma: no cover

    parser = argparse.ArgumentParser(description="Loads vertexes and edges into postgres")
    parser.add_argument('--dry', dest='dry', action='store_true', default=False, help="echo commands, don't execute [False]")
    parser.add_argument('--drop', dest='drop', action='store_true', default=False, help="drop the tables first [False]")
    parser.add_argument('--skip_index', dest='index', action='store_false', default=True, help="index after loading [True]")
    config_path = "{}/config.yml".format(os.path.dirname(os.path.realpath(__file__)))
    parser.add_argument('--config', dest='config', default=config_path, help="config path {}".format(config_path))
    args = parser.parse_args()
    print(vars(args))
    main(**vars(args))
