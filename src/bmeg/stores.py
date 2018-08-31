
"""
store and retrieve anything with a gid()
"""
from bmeg.vertex import Allele, AlleleAnnotations
from bmeg.ioutils import reader
import dataclasses
import logging
import sqlite3
import os.path

import json


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
            "insert into data values(?, ?, ?)",
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

    def get(self, gid):
        """ xform dict to Allele"""
        return Allele.from_dict(super(MyVariantSqlitestore, self).get(gid))


def new_store(name, **kwargs):
    """ create store based on names"""
    if name == 'memory':
        return Memorystore()
    if name == 'myvariantinfo-memory':
        return MyVariantTestMemorystore(**kwargs)
    if name == 'myvariantinfo-sqlite':
        return MyVariantSqlitestore(**kwargs)
    assert False, 'no store named {}'.format(name)
