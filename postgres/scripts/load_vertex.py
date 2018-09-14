"""
loads json vertex files into postgres
"""

import dataset
import ujson
import logging

FILES = """
outputs/gtex/gtex.Project.Vertex.json
outputs/gtex/gtex.Biosample.Vertex.json
""".strip().split()

BATCH_SIZE = 1000
logging.getLogger().setLevel(logging.INFO)


def construct_pg_url(postgres_user='postgres', postgres_password='password', postgres_host='database', postgres_port='5432', postgres_database='test'):
    PG_URL = "postgresql://" + postgres_user + ":" + postgres_password + '@' + postgres_host + ':' + postgres_port + '/' + postgres_database
    return PG_URL


pgconn = dataset.Database(url=construct_pg_url())


def rows():
    t = c = 0
    for f in FILES:
        logging.info('reading {}'.format(f))
        with open(f) as ins:
            for line in ins:
                c += 1
                t += 1
                obj = ujson.loads(line)
                del obj['_id']
                yield obj
                if c % BATCH_SIZE == 0:
                    c = 0
                    logging.info('loaded {}'.format(t))


vertex = pgconn['vertex']
vertex.insert_many(rows())

logging.info('There are {} vertexes'.format(vertex.count()))
