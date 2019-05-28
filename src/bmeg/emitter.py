import atexit
import json
import os
import sys
import typing
import dataclasses
from datetime import datetime
import gzip
from collections import Iterable

from bmeg import ClassInstance
from bmeg.gid import GID
from bmeg.utils import enforce_types, ensure_directory


class Rate:
    def __init__(self):
        self.i = 0
        self.start = None
        self.first = None

    def close(self):
        if self.i == 0:
            return

        self.log()
        dt = datetime.now() - self.first
        m = "\ntotal: {0:,} in {1:,d} seconds".format(self.i,
                                                      int(dt.total_seconds()))
        print(m, file=sys.stderr)

    def log(self):
        if self.i == 0:
            return

        dt = datetime.now() - self.start
        self.start = datetime.now()
        rate = 1000 / dt.total_seconds()
        m = "rate: {0:,} emitted ({1:,d}/sec)".format(self.i, int(rate))
        print("\r" + m, end='', file=sys.stderr)

    def tick(self):
        if self.start is None:
            self.start = datetime.now()
            self.first = self.start

        self.i += 1

        if self.i % 1000 == 0:
            self.log()

class DebugEmitter:
    def __init__(self, **kwargs):
        self.generator = Generator(**kwargs)

    def close(self):
        self.generator.close()

    def write_edge(self, e):
        print(json.dumps(e, indent=True))

    def write_vertex(self, v):
        print(json.dumps(v, indent=True))

    def emit(self, obj):
        o = self.generator.to_dict(obj)
        print(json.dumps(o, indent=True))

    def emit_edge(self, obj, from_gid: GID, to_gid: GID):
        self.generator.emit_edge(
            self,
            obj, from_gid, to_gid
        )

    def emit_vertex(self, obj):
        self.generator.emit_vertex(
            self,
            obj
        )

    def emit_link(self, label, to_gid, from_gid):
        self.write_edge({"label" : label, "to" : to_gid, "from" : from_gid})


class JSONEmitter:
    def __init__(self, directory, prefix=None, **kwargs):
        self.vertex_handles = FileHandler(directory, prefix, "Vertex.json")
        self.edge_handles = FileHandler(directory, prefix, "Edge.json")
        self.obj_handles = FileHandler(directory, prefix, "json")
        self.generator = Generator(**kwargs)

    def close(self):
        self.vertex_handles.close()
        self.edge_handles.close()
        self.generator.close()

    def write_edge(self, obj):
        fh = self.edge_handles[obj]
        if self.edge_handles.compresslevel > 0:
            fh.write(json.dumps(obj).encode())
            fh.write(os.linesep.encode())
        else:
            fh.write(json.dumps(obj))
            fh.write(os.linesep)

    def write_vertex(self, obj):
        fh = self.vertex_handles[obj]
        if self.vertex_handles.compresslevel > 0:
            fh.write(json.dumps(obj).encode())
            fh.write(os.linesep.encode())
        else:
            fh.write(json.dumps(obj))
            fh.write(os.linesep)

    def write_obj(self, obj):
        o = self.generator.to_dict(obj)
        fh = self.obj_handles[obj]
        if self.vertex_handles.compresslevel > 0:
            fh.write(json.dumps(o).encode())
            fh.write(os.linesep.encode())
        else:
            fh.write(json.dumps(o))
            fh.write(os.linesep)

    def emit(self, obj):
        self.write_obj(obj)

    def emit_edge(self, obj, from_gid: GID, to_gid: GID):
        d = self.generator.emit_edge(self, obj, from_gid, to_gid)

    def emit_vertex(self, obj):
        d = self.generator.emit_vertex(self, obj)

    def emit_link(self, label, to_gid, from_gid):
        self.write_edge({"label" : label, "to" : to_gid, "from" : from_gid})

def make_gid(label, from_gid: GID, to_gid: GID):
    return "(%s)--%s->(%s)" % (from_gid, label, to_gid)

class Generator:
    """
    Generator is an internal helper that contains code shared by all
    emitters, such as validation checks, data cleanup, etc.
    """

    def __init__(self, preserve_null=False):
        self.preserve_null = preserve_null
        self.rate = Rate()

    def close(self):
        self.rate.close()

    def _get_data(self, obj):
        # this util recurses and unravels embedded dataclasses
        # see https://docs.python.org/3/library/dataclasses.html#dataclasses.asdict
        data = dataclasses.asdict(obj)
        # delete null values
        if not self.preserve_null:
            remove = [k for k in data if data[k] is None]
            for k in remove:
                del data[k]

        return data

    @enforce_types
    def to_dict(self, obj: ClassInstance):
        err = obj._validate()
        if len(err) > 0:
            #TODO: emit errors to log
            for i in err:
                print("Error: %s" % (i))
            #return None
        return obj._data

    @enforce_types
    def emit_edge(self, wrt, obj: ClassInstance, from_gid: str, to_gid: str):
        label = obj._classSchema.label
        gid = make_gid(label, from_gid, to_gid)
        data = {}
        for k, val in obj._data.items():
            link = obj._classSchema.getLink(k)
            if link is None:
                data[k] = val

        dumped = {
            "_id": gid,
            "gid": gid,
            "label": label,
            "from": from_gid,
            "to": to_gid,
            "data": data
        }
        self.rate.tick()
        wrt.write_edge( dumped )

    @enforce_types
    def emit_vertex(self, wrt, obj: ClassInstance):
        err = obj._validate()
        if len(err) > 0:
            #TODO: emit errors to log
            for i in err:
                print("Error: %s" % (i))
            #return None

        gid = None
        label = obj._classSchema.label
        for k, prop in obj._classSchema.properties().items():
            # If the value represents the 'node_id' of the data
            if prop.get("systemAlias", "") == "node_id":
                gid = obj._get(k)
        if gid is None:
            print("%s node_id is not set" % (obj._classSchema.name))
            return None
        data = {}
        for k, val in obj._data.items():
            link = obj._classSchema.getLink(k)
            if link is not None:
                if not isinstance(val, Iterable) or isinstance(val, str):
                    val = [val]
                for dst in val:
                    self.emit_link(wrt, link['label'], from_gid=gid, to_gid=dst)
                    if 'backref' in link:
                        self.emit_link(wrt, link["backref"], to_gid=gid, from_gid=dst)
            #elif obj._classSchema.prop(k).systemAlias == "node_id":
            #    pass
            else:
                data[k] = val
        wrt.write_vertex( {"gid":gid, "label":label, "data" : data} )
        self.rate.tick()

    # a link if an edge with no properties
    def emit_link(self, wrt, label, from_gid, to_gid ):
        wrt.write_edge({"label":label, "from":from_gid, "to":to_gid})


class FileHandler:
    """
    FileHandler helps manage a set of file handles, indexed by a key.
    This is used by emitters to write to a set of files, such as
    Sample.Vertex.json, Case.Vertex.json, etc.

    This is an internal helper.
    """
    def __init__(self, directory, prefix, extension, mode="w", compresslevel=1):
        self.prefix = prefix
        self.directory = directory
        ensure_directory("outputs", self.directory)
        self.outdir = os.path.join("outputs", self.directory)
        self.extension = extension
        self.mode = mode
        self.handles = {}
        self.compresslevel = compresslevel
        atexit.register(self.close)

    def __getitem__(self, obj):
        if isinstance(obj, ClassInstance):
            label = obj._classSchema.name
        else:
            label = obj['label'] #.__class__.__name__

        fname = "%s.%s" % (label, self.extension)
        if self.prefix is not None:
            fname = self.prefix + "." + fname
        fname = os.path.join(self.outdir, fname)

        if fname in self.handles:
            return self.handles[fname]
        else:
            if self.compresslevel:
                self.mode = "wb"
                fh = gzip.GzipFile(
                    filename='',
                    compresslevel=self.compresslevel,
                    fileobj=open(fname + '.gz', mode=self.mode),
                    mtime=0
                )
            else:
                fh = open(fname, self.mode)
            self.handles[fname] = fh
            return fh

    def close(self):
        for fh in self.handles.values():
            fh.close()


# shorthand aliases for emitter names
EMITTER_NAME_MAP = {
    "json": JSONEmitter,
    "debug": DebugEmitter,
}


def new_emitter(name="json", **kwargs):
    """ construct an emitter """
    cls = EMITTER_NAME_MAP[name]
    return cls(**kwargs)
