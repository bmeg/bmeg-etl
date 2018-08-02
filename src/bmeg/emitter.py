import atexit
import json
import os
import sys
import typing

from datetime import datetime

from bmeg.edge import Edge
from bmeg.gid import GID
from bmeg.utils import enforce_types
from bmeg.vertex import Vertex


class DebugEmitter:
    def __init__(self, **kwargs):
        self.emitter = BaseEmitter(**kwargs)

    def close(self):
        self.emitter.close()

    def emit_edge(self, obj: Edge, from_gid: GID, to_gid: GID):
        d = self.emitter.emit_edge(obj)
        print(json.dumps(d, indent=True))

    def emit_vertex(self, obj: Vertex):
        d = self.emitter.emit_vertex(obj)
        print(json.dumps(d, indent=True))


class JSONEmitter:
    def __init__(self, prefix, **kwargs):
        self.handles = FileHandler(prefix, "json")
        self.emitter = BaseEmitter(**kwargs)

    def close(self):
        self.handles.close()
        self.emitter.close()

    def emit_edge(self, obj: Edge, from_gid: GID, to_gid: GID):
        d = self.emitter.emit_edge(obj)
        fh = self.handles[obj]
        json.dump(d, fh)
        fh.write(os.linesep)

    def emit_vertex(self, obj: Vertex):
        d = self.emitter.emit_vertex(obj)
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


class BaseEmitter:
    """
    BaseEmitter is an internal helper that contains code shared by all
    emitters, such as validation checks, data cleanup, etc.
    """

    def __init__(self, preserve_null=False):
        self.preserve_null = preserve_null
        self.rate = Rate()

    def close(self):
        self.rate.close()

    def _get_data(self, obj: typing.Union[Edge, Vertex]):
        data = dict(obj.__dict__)

        # delete null values
        if not self.preserve_null:
            remove = [k for k in data if data[k] is None]
            for k in remove:
                del data[k]

        return data

    @enforce_types
    def emit_edge(self, obj: Edge, from_gid: GID, to_gid: GID):
        dumped = {
            "gid": obj.make_gid(from_gid, to_gid),
            "label": obj.label(),
            "from": from_gid,
            "to": to_gid,
            "data": self._get_data(obj)
        }

        self.rate.tick()
        return dumped

    @enforce_types
    def emit_vertex(self, obj: Vertex):
        dumped = {
            "gid": obj.gid(),
            "label": obj.label(),
            "data": self._get_data(obj)
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
