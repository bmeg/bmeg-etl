
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



class Sqlitestore:

    def __init__(self, path):
        self.conn = sqlite3.connect(path)
        # optimize db calls
        self.conn.execute("PRAGMA synchronous = OFF;")
        self.conn.execute("PRAGMA journal_mode = OFF;")
        # create table
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

    def all_ids(self):
        cur = self.conn.cursor()
        for t in cur.execute('SELECT gid FROM data ;'):
            yield t[0]


class AlleleSqlitestore(Sqlitestore):

    def get(self, gid):
        """ xform dict to Allele"""
        return Allele.from_dict(super(AlleleSqlitestore, self).get(gid))

    def all(self):
        """ xform dict to Allele"""
        for allele in super(AlleleSqlitestore, self).all():
            yield Allele.from_dict(allele)

    def load_many(self, batch):
        """ load a batch of tuples into table """
        logging.info('starting insert')
        self.conn.executemany("insert or replace into data(gid, clazz, json) values(?, ?, ?)", batch)
        logging.info('done insert')


class KeyValuestore:
    """ store an id an a json serializable object"""
    def __init__(self, path):
        self.conn = sqlite3.connect(path)
        # optimize db calls
        self.conn.execute("PRAGMA synchronous = OFF;")
        self.conn.execute("PRAGMA journal_mode = OFF;")
        # create table
        cur = self.conn.cursor()
        cur.execute("CREATE TABLE IF NOT EXISTS data (id text, json text);")
        self.conn.commit()

    def index(self):
        cur = self.conn.cursor()
        cur.execute("CREATE UNIQUE INDEX IF NOT EXISTS idx_data ON data(id);")
        self.conn.commit()

    def get(self, id):
        result = None
        cur = self.conn.cursor()
        cur.execute("select json from data where id=?", (id,))
        t = cur.fetchone()
        if t:
            result = json.loads(t[0])
        cur.close()
        return result

    def put(self, id, obj):
        cur = self.conn.cursor()
        cur.execute(
            "insert or replace into data values(?, ?)",
            (id, json.dumps(obj, separators=(',', ':')))
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
        for t in cur.execute('SELECT json FROM data ;'):
            yield json.loads(t[0])

    def all_ids(self):
        cur = self.conn.cursor()
        for t in cur.execute('SELECT id FROM data ;'):
            yield t[0]

    def load_many(self, batch):
        """ load a batch of tuples into table """
        logging.info('starting insert')
        self.conn.executemany("insert or replace into data(id, json) values(?, ?)", batch)
        logging.info('done insert')


def new_store(name, **kwargs):
    """ create store based on names"""
    if name == 'memory':
        return Memorystore()
    if name == 'sqlite':
        return Sqlitestore(**kwargs)
    if name == 'allele-sqlite':
        return AlleleSqlitestore(**kwargs)
    if name == 'compound':
        return KeyValuestore(**kwargs)
    assert False, 'no store named {}'.format(name)
