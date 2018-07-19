from bmeg.models.vertex_models import Vertex
from bmeg.models.edge_models import Edge


class Emitter:
    def __init__(self, prefix):
        self.prefix = prefix
        self.handles = {}

    def close(self):
        for fname, fh in self.handles.items():
            fh.close()

    def emit(self, obj):
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

        fh.write(obj.dump())
        return
