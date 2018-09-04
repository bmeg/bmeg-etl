
"""
store and retrieve anything with a gid()
"""
from bmeg.vertex import Allele, AlleleAnnotations
from bmeg.ioutils import reader
import dataclasses
import logging
import sqlite3
import os.path
import hashlib
import json
from threading import Thread
from queue import Queue
import queue


class Memorystore:

    def __init__(self):
        self.key_val = {}

    def get(self, gid):
        return self.key_val.get(gid, None)

    def put(self, obj):
        self.key_val[obj.gid()] = obj

    def all(self):
        for k in self.key_val.keys():
            yield self.key_val[k]

    def size(self):
        return len(self.key_val.keys())


class MyVariantTestMemorystore(Memorystore):

    def __init__(self, myvariantinfo_path):
        super(MyVariantTestMemorystore, self).__init__()
        self.load(myvariantinfo_path)

    def load(self, myvariantinfo_path):
        with reader(myvariantinfo_path) as ins:
            for line in ins:
                myvariant = json.loads(line)
                aa = AlleleAnnotations()
                aa.myvariantinfo = myvariant
                allele = Allele(
                    genome='GRCh37',
                    chromosome=myvariant['chrom'],
                    start=myvariant['hg19']['start'],
                    end=myvariant['hg19']['end'],
                    reference_bases=myvariant['vcf']['ref'],
                    alternate_bases=myvariant['vcf']['alt'],
                    annotations=aa
                )
                self.put(allele)


class Sqlitestore:

    def __init__(self, path='sqllite.db'):
        self.conn = sqlite3.connect(path)
        cur = self.conn.cursor()
        cur.execute("CREATE TABLE IF NOT EXISTS data (gid text, clazz text, json text);")
        self.conn.commit()


    def index(self):
        cur = self.conn.cursor()
        cur.execute("CREATE UNIQUE INDEX IF NOT EXISTS idx_data ON data(gid);")
        self.conn.commit()

    def get(self, gid):
        result = None
        cur = self.conn.cursor()
        cur.execute("select * from data where gid=?", (gid,))
        t = cur.fetchone()
        if t:
            result = json.loads(t[2])
        cur.close()
        return result

    def put(self, obj):
        cur = self.conn.cursor()
        cur.execute(
            "insert or replace into data values(?, ?, ?)",
            (obj.gid(), obj.__class__.__name__, json.dumps(dataclasses.asdict(obj), separators=(',', ':')))
        )
        self.conn.commit()
        cur.close()

    def size(self):
        result = 0
        try:
            cur = self.conn.cursor()
            result = cur.execute("select MAX(_ROWID_) from data LIMIT 1;").fetchone()[0]
            cur.close()
        except Exception:
            pass
        if not result:
            result = 0
        return result

    def all(self):
        cur = self.conn.cursor()
        for t in cur.execute('SELECT * FROM data ;'):
            yield json.loads(t[2])


class AlleleSqlitestore(Sqlitestore):

    def get(self, gid):
        """ xform dict to Allele"""
        return Allele.from_dict(super(AlleleSqlitestore, self).get(gid))

    def all(self, gid):
        """ xform dict to Allele"""
        for allele in super(AlleleSqlitestore, self).all():
            yield Allele.from_dict(allele)


class MyVariantSqlitestore(Sqlitestore):
    """ create sqlite.db in same directory as myvariantinfo_path """
    def __init__(self, myvariantinfo_path):
        path = '{}/{}'.format(os.path.dirname(myvariantinfo_path), 'sqlite.db')
        super(MyVariantSqlitestore, self).__init__(path=path)
        size = self.size()
        if size == 0:
            logging.warning('loading. this may take a while {}'.format(path))
            self.load(myvariantinfo_path)
        else:
            logging.info('proceeding. there are {} rows in table'.format(size))

    def load(self, myvariantinfo_path):
        """ take myvariant into json, xform to Allele"""
        def make_gid(cls, genome, chromosome, start, end, reference_bases,
                     alternate_bases):
            vid = "%s:%s:%d:%d:%s:%s" % (genome, chromosome,
                                         start, end, reference_bases,
                                         alternate_bases)
            vid = vid.encode('utf-8')
            vidhash = hashlib.sha1()
            vidhash.update(vid)
            vidhash = vidhash.hexdigest()
            return "%s:%s" % (cls, vidhash)



        def reader_worker():
            t = c = 0
            batch = []
            batch_size = 20000
            with reader(myvariantinfo_path) as ins:
                for line in ins:
                    myvariant = json.loads(line)
                    if 'hg19' not in myvariant:
                        continue
                    if 'vcf' not in myvariant:
                        continue
                    # aa = AlleleAnnotations()
                    # aa.myvariantinfo = myvariant
                    # allele = Allele(
                        # genome='GRCh37',
                        # chromosome=myvariant['chrom'],
                        # start=myvariant['hg19']['start'],
                        # end=myvariant['hg19']['end'],
                        # reference_bases=myvariant['vcf']['ref'],
                        # alternate_bases=myvariant['vcf']['alt'],
                    #     annotations=aa
                    # )
                    # batch.append(
                    #     (
                    #         allele.gid(),
                    #         allele.__class__.__name__,
                    #         json.dumps(dataclasses.asdict(allele), separators=(',', ':'))
                    #     )
                    # )
                    gid = make_gid(
                        cls='Allele',
                        genome='GRCh37',
                        chromosome=myvariant['chrom'],
                        start=myvariant['hg19']['start'],
                        end=myvariant['hg19']['end'],
                        reference_bases=myvariant['vcf']['ref'],
                        alternate_bases=myvariant['vcf']['alt'],
                    )
                    allele = {
                        'genome': 'GRCh37',
                        'chromosome': myvariant['chrom'],
                        'start': myvariant['hg19']['start'],
                        'end': myvariant['hg19']['end'],
                        'reference_bases': myvariant['vcf']['ref'],
                        'alternate_bases': myvariant['vcf']['alt'],
                        'annotations': {'myvariantinfo': myvariant}
                    }
                    batch.append(
                        (
                            gid,
                            'Allele',
                            json.dumps(allele, separators=(',', ':'))
                        )
                    )

                    c += 1
                    t += 1
                    if c % batch_size == 0:
                        c = 0
                        logging.info('loaded {}'.format(t))
                        pending_q.put(batch)
                        batch = []
                pending_q.put('DONE')
                logging.info('loaded {}'.format(t))

        # q to hold messages
        pending_q = Queue(maxsize=5)
        logging.debug('created pending_q')
        # optimize db calls
        self.conn.execute("PRAGMA synchronous = OFF;")
        self.conn.execute("PRAGMA journal_mode = OFF;")

        # start worker
        worker = Thread(target=reader_worker, args=(), name='reader_worker')
        worker.setDaemon(True)
        worker.start()

        while True:
            logging.info('waitin on q')
            batch = pending_q.get(block=True, timeout=20)
            if not batch:
                contine
            if batch == 'DONE':
                break
            logging.info('starting insert')
            self.conn.executemany("insert into data(gid, clazz, json) values(?, ?, ?)", batch)
            logging.info('done insert')
        logging.info('creating index')
        self.index()
        logging.info('done creating index')




    def get(self, gid):
        """ xform dict to Allele"""
        return Allele.from_dict(super(MyVariantSqlitestore, self).get(gid))

    def all(self, gid):
        """ xform dict to Allele"""
        for allele in super(MyVariantSqlitestore, self).all():
            yield Allele.from_dict(allele)


def new_store(name, **kwargs):
    """ create store based on names"""
    if name == 'memory':
        return Memorystore()
    if name == 'sqlite':
        return Sqlitestore()
    if name == 'myvariantinfo-memory':
        return MyVariantTestMemorystore(**kwargs)
    if name == 'myvariantinfo-sqlite':
        return MyVariantSqlitestore(**kwargs)
    if name == 'allele-sqlite':
        return AlleleSqlitestore(**kwargs)
    assert False, 'no store named {}'.format(name)
