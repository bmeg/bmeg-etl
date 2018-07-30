import atexit
from datetime import datetime
import json
import os
import sys

from bmeg.models.vertex_models import Vertex
from bmeg.models.edge_models import Edge


class DebugEmitter:
    def __init__(self, **kwargs):
        self.emitter = BaseEmitter(**kwargs)

    def close(self):
        self.emitter.close()

    def emit(self, obj):
        d = self.emitter.emit(obj)
        print(json.dumps(d, indent=True))


class JSONEmitter:
    def __init__(self, prefix, **kwargs):
        self.handles = FileHandler(prefix, "json")
        self.emitter = BaseEmitter(**kwargs)

    def close(self):
        self.handles.close()
        self.emitter.close()

    def emit(self, obj):
        d = self.emitter.emit(obj)
        fh = self.handles[obj]
        json.dump(d, fh)
        fh.write(os.linesep)


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
        m = "\ntotal: {0:,} in {1:,d} seconds".format(self.i, int(dt.total_seconds()))
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


class BaseEmitter:
    """
    BaseEmitter is an internal helper that contains code shared by all emitters,
    such as validation checks, data cleanup, etc.
    """

    def __init__(self, preserve_null=False):
        self.preserve_null = preserve_null
        self.rate = Rate()

    def close(self):
        self.rate.close()

    def emit(self, obj):

        if not obj.gid:
            raise ValueError("gid is empty")

        if not isinstance(obj, Vertex) and not isinstance(obj, Edge):
            raise TypeError("emit accepts objects of the Vertex or Edge type")

        label = obj.__class__.__name__
        data = dict(obj.__dict__)

        if "gid" in data:
            del data["gid"]

        if "from_gid" in data:
            del data["from_gid"]

        if "to_gid" in data:
            del data["to_gid"]

        # delete null values
        if not self.preserve_null:
            remove = [k for k in data if data[k] is None]
            for k in remove:
                del data[k]

        if isinstance(obj, Vertex):
            dumped = {
                "gid": obj.gid,
                "label": label,
                "data": data
            }

        elif isinstance(obj, Edge):
            suffix = "Edge"
            dumped = {
                "gid": obj.gid,
                "label": label,
                "from": obj.from_gid,
                "to": obj.to_gid,
                "data": data
            }

        self.rate.tick()
        return dumped


class FileHandler:
    """
    FileHandler helps manage a set of file handles, indexed by a key.
    This is used by emitters to write to a set of files, such as
    Biosample.Vertex.json, Individual.Vertex.json, etc.

    This is an internal helper.
    """
    def __init__(self, prefix, extension, mode="w"):
        self.prefix = prefix
        self.extension = extension
        self.mode = mode
        self.handles = {}
        atexit.register(self.close)

    def __getitem__(self, obj):
        label = obj.__class__.__name__

        if isinstance(obj, Vertex):
            suffix = "Vertex"
        elif isinstance(obj, Edge):
            suffix = "Edge"
        else:
            suffix = "Unknown"

        fname = "%s.%s.%s.%s" % (self.prefix, label, suffix, self.extension)

        if fname in self.handles:
            return self.handles[fname]
        else:
            fh = open(fname, self.mode)
            self.handles[fname] = fh
            return fh
        
    def close(self):
        for fh in self.handles.values():
            fh.close()


