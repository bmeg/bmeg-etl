"""
loads json vertex files into postgres table edge(gid, label, from, to, data)
"""

import dataset
import ujson
import logging
import yaml
import types
from jinja2 import Template
import os
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

pgconn = dataset.Database(url=construct_pg_url(**config.postgres))


def execute(pgconn, commands):
    """ exec statement on connection and commit """
    pgconn.begin()
    for command in commands:
        logging.info(command)
        pgconn.query(command)
    pgconn.commit()


def rows(files, batch_size=1000):
    """
    generator: read in all rows from all files,
    remove keys_to_delete before yielding
    log message every batch_size
    """
    for f in files:
        logging.info('reading {}'.format(f))
        t = c = 0
        with reader(f) as ins:
            for line in ins:
                c += 1
                t += 1
                obj = ujson.loads(line)
                # if edge, copy from and to underneath data
                obj['data']['gid'] = obj['gid']
                if 'from' in obj:
                    obj['data']['from'] = obj['from']
                if 'to' in obj:
                    obj['data']['to'] = obj['to']
                yield obj
                if c % batch_size == 0:
                    c = 0
                    logging.info('loaded {}'.format(t))


def vertex_ddl(files, vertex_template):
    """
    generator: for every vertex file, create cypher vertex
    """
    t = Template(vertex_template)
    for f in files:
        logging.info('reading {}'.format(f))
        # read first line to get label
        label = None
        with reader(f) as ins:
            for line in ins:
                label = ujson.loads(line)['label']
                break
        yield t.render(vertex=label)


def edge_ddl(files, edge_template):
    """
    generator: for every edge file, create cypher edge
    """
    t = Template(edge_template)
    for f in files:
        logging.info('reading {}'.format(f))
        # read first line to get label
        label = from_label = to_label = None
        with reader(f) as ins:
            for line in ins:
                label = ujson.loads(line)['label']
                from_label = ujson.loads(line)['from'].split(':')[0]
                to_label = ujson.loads(line)['to'].split(':')[0]
                break
        yield t.render(edge=label, from_label=from_label, to_label=to_label)


logging.info('(re) creating tables')
execute(pgconn, [config.ddl])

# connect to table, load it and log count

vextexes = pgconn['bmeg_vertex']
logging.info('creating vertexes table for loading')
vextexes.insert_many(rows(config.vertex_files), chunk_size=100)
logging.info('There are {} vextexes'.format(vextexes.count()))

logging.info('creating graph vertexes')
execute(pgconn, [ddl for ddl in vertex_ddl(config.vertex_files, config.vertex_template)])

edges = pgconn['bmeg_edge']
logging.info('creating edge table for loading')
edges.insert_many(rows(config.edge_files), chunk_size=100)
logging.info('There are {} edges'.format(edges.count()))

logging.info('creating graph edges')
execute(pgconn, [ddl for ddl in edge_ddl(config.edge_files, config.edge_template)])

# logging.info('creating indexes')
# execute(pgconn, [config.indexes])
# logging.info('indexes created')
