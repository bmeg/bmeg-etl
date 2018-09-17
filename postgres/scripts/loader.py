"""
loads json vertex files into postgres table edge(gid, label, from, to, data)
"""

import dataset
import ujson
import logging
import yaml


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
# host: database
# user: postgres
# passwd: password
# db: test
def construct_pg_url(user='postgres', password='password', host='database', port='5432', database='test'):
    PG_URL = "postgresql://" + user + ":" + password + '@' + host + ':' + port + '/' + database
    return PG_URL


with open("config.yml", 'r') as stream:
    config = yaml.load(stream)

pgconn = dataset.Database(url=construct_pg_url())


def execute(pgconn, commands):
    """ exec statement on connection and commit """
    pgconn.begin()
    for command in commands:
        pgconn.query(command)
    pgconn.commit()


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
                yield obj
                if c % batch_size == 0:
                    c = 0
                    logging.info('loaded {}'.format(t))


logging.info('(re) creating tables')
execute(pgconn, [DROP_TABLES, CREATE_TABLES])

# connect to table, load it and log count
edges = pgconn['edge']
logging.info('inserting edges')
edges.insert_many(rows(EDGE_FILES))
logging.info('There are {} edges'.format(edges.count()))

vextexes = pgconn['vertex']
logging.info('inserting vertexes')
vextexes.insert_many(rows(VERTEX_FILES))
logging.info('There are {} vextexes'.format(vextexes.count()))

logging.info('creating indexes')
execute(pgconn, [CREATE_INDEXES])
