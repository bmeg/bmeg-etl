
"""
store and retrieve anything with an id
"""
import logging
import sqlite3
import ujson
import uuid

from bmeg import ClassInstance


class KeyValueMemoryStore:
    """ store an id and a json serializable object in sqllite"""

    def __init__(self):
        self.key_val = {}

    def backend(self):
        """ return the implementation dependent backend"""
        return self.key_val

    def get(self, id):
        """ return single obj, none if no match """
        return self.key_val.get(id, None)

    def put(self, id, obj):
        """ save obj """
        self.key_val[id] = obj

    def delete(self, id):
        """ remove obj """
        del self.key_val[id]

    def all(self):
        """ yield all """
        for k in self.key_val.keys():
            yield self.key_val[k]

    def all_ids(self):
        """ yield all ids """
        for k in self.key_val.keys():
            yield k

    def size(self):
        """ the number of objects stored """
        return len(self.key_val.keys())

    def load_many(self, batch):
        """ load a batch of tuples into table. each tuple is (str, obj) """
        for b in batch:
            self.put(id=b[0], obj=b[1])


class KeyValueStore:
    """ store an id and a json serializable object in sqllite"""
    def __init__(self, path=None, index=False):
        if not path:
            path = '/tmp/{}.db'.format(uuid.uuid4())
        self.conn = sqlite3.connect(path)
        # optimize db calls
        self.conn.execute("PRAGMA synchronous = OFF;")
        self.conn.execute("PRAGMA journal_mode = OFF;")
        # create table
        cur = self.conn.cursor()
        cur.execute("CREATE TABLE IF NOT EXISTS data (id text, json text);")
        self.conn.commit()
        if index:
            self.index()

    def backend(self):
        """ return the implementation dependent backend"""
        return self.conn

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
            result = ujson.loads(t[0])
        cur.close()
        return result

    def put(self, id, obj):
        if isinstance(obj, ClassInstance):
            obj = obj.props()
        cur = self.conn.cursor()
        cur.execute(
            "insert or replace into data values(?, ?)",
            (id, ujson.dumps(obj))
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
            logging.info('could not read max(_rowid_)', exc_info=True)
        if not result:
            logging.info('max(_rowid_) returned null')
            result = 0
        return result

    def all(self):
        cur = self.conn.cursor()
        for t in cur.execute('SELECT json FROM data ;'):
            yield ujson.loads(t[0])

    def all_ids(self):
        cur = self.conn.cursor()
        for t in cur.execute('SELECT id FROM data ;'):
            yield t[0]

    def load_many(self, batch):
        """ load a batch of tuples into table. each tuple is (str, json-serializable) """
        # serialize batch
        pbatch = []
        for idx, i in enumerate(batch):
            if isinstance(i[1], ClassInstance):
                pbatch.append((i[0], i[1].props()))
            else:
                pbatch.append((i[0], i[1]))
        batch = [(i[0], ujson.dumps(i[1])) for i in pbatch]
        logging.debug('starting insert')
        self.conn.executemany("insert or replace into data(id, json) values(?, ?)", batch)
        logging.debug('done insert')


class DataClassStore(KeyValueStore):
    """ stores a dataclass, uses gid as id, loads dataobject """
    def __init__(self, clazz, **kwargs):
        super().__init__(**kwargs)
        self.clazz = clazz

    def get(self, gid):
        """ xform dict to dataclass"""
        obj = super().get(gid)
        return self.clazz(**obj)

    def put(self, obj):
        """ get gid from dataclass"""
        return super().put(obj.gid(), obj)

    def all(self):
        """ xform dict to dataclass"""
        for obj in super().all():
            yield self.clazz(**obj)

    def load_many(self, batch):
        """ load a batch of tuples into table. each tuple is (str, object) """
        batch = [(b.gid(), b) for b in batch]
        super().load_many(batch)


class DataClassMemoryStore(KeyValueMemoryStore):
    """ stores a dataclass, uses gid as id, loads dataobject """
    def __init__(self, clazz):
        super().__init__()
        self.clazz = clazz

    def put(self, obj):
        """ get gid from dataclass"""
        return super().put(obj.gid(), obj)

    def load_many(self, batch):
        """ load a batch of tuples into table. each tuple is (str, object) """
        for b in batch:
            super().put(id=b.gid(), obj=b)


def new_store(name, **kwargs):
    """ create store based on names"""
    if name == 'memory':
        return KeyValueMemoryStore()
    if name == 'key-val':
        return KeyValueStore(**kwargs)
    if name == 'dataclass':
        return DataClassStore(**kwargs)
    if name == 'dataclass-memory':
        return DataClassMemoryStore(**kwargs)
    raise NotImplementedError('no store named {}'.format(name))
