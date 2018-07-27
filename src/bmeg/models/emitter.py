import json
import os

from bmeg.models.vertex_models import Vertex
from bmeg.models.edge_models import Edge


class Emitter:
    def __init__(self, prefix, preserve_null=False):
        self.prefix = prefix
        self.handles = {}
        self.preserve_null = preserve_null

    def close(self):
        for fname, fh in self.handles.items():
            fh.close()

    def emit(self, obj):

        if not obj.gid:
            raise ValueError("gid is empty")

        label = obj.__class__.__name__

        if isinstance(obj, Vertex):
            suffix = "Vertex"
        elif isinstance(obj, Edge):
            suffix = "Edge"
        else:
            raise TypeError("emit accepts objects of the Vertex or Edge type")

        fname = "%s.%s.%s.json" % (self.prefix, label, suffix)

        if fname in self.handles:
            fh = self.handles[fname]
        else:
            fh = open(fname, "w")
            self.handles[fname] = fh

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
            dumped = json.dumps({
                "gid": self.gid,
                "label": label,
                "data": data
            })

        elif isinstance(obj, Edge):
            suffix = "Edge"
            dumped = json.dumps({
                "gid": obj.gid,
                "label": label,
                "from": obj.from_gid,
                "to": obj.to_gid,
                "data": data
            })

        fh.write(dumped)
        fh.write(os.linesep)

        return
