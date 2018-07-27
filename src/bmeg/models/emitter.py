import json
import os

from bmeg.models.vertex_models import Vertex
from bmeg.models.edge_models import Edge


class JSONEmitter:
    def __init__(self, prefix, **kwargs):
        self.handles = filehandler(prefix)
        self.emitter = emitter(**kwargs)

    def close(self):
        self.handles.close()

    def emit(self, obj):
        d = self.emitter.emit(obj)
        fh = self.handles[obj]
        json.dump(d, fh)
        fh.write(os.linesep)


class emitter:
    """
    emitter is an internal helper that contains code shared by all emitters,
    such as validation checks, data cleanup, etc.
    """

    def __init__(self, preserve_null=False):
        self.preserve_null = preserve_null

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
        if not self.preserve_nulls:
            remove = [k for k in data if data[k] is None]
            for k in remove:
                del data[k]

        if isinstance(obj, Vertex):
            dumped = {
                "gid": self.gid,
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

        return dumped


class filehandler:
    """
    filehandler helps manage a set of file handles, indexed by a key.
    This is used by emitters to write to a set of files, such as
    Biosample.Vertex.json, Individual.Vertex.json, etc.

    This is an internal helper.
    """
    def __init__(self, prefix):
        self.prefix = prefix
        self.handles = {}

    def __get__(self, obj):
        label = obj.__class__.__name__

        if isinstance(obj, Vertex):
            suffix = "Vertex"
        elif isinstance(obj, Edge):
            suffix = "Edge"
        else:
            suffix = "Unknown"

        fname = "%s.%s.%s.json" % (self.prefix, label, suffix)

        if fname in self.handles:
            return self.handles[fname]
        else:
            fh = open(fname, "w")
            self.handles[fname] = fh
            return fh
        
    def close(self):
        for fh in self.handles.values():
            fh.close()


