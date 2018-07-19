import json

from dataclasses import dataclass, field

from bmeg.models.gid_models import GID
from bmeg.models.utils import set_gid


@dataclass(frozen=True)
class Edge:
    gid: GID = field(init=False)

    def dump(self):
        return json.dumps({
            "gid": self.gid,
            "label": self.__class__.__name__,
            "from": self.from_gid,
            "to": self.to_gid,
            "data": self.data
        })


@dataclass(frozen=True)
class GenericEdge(Edge):
    label: str
    from_gid: GID
    to_gid: GID

    def __post_init__(self):
        set_gid(self, Edge.make_gid(self.label, self.from_gid, self.to_gid))

    @classmethod
    def make_gid(cls, label, from_gid, to_gid):
        return "(%s)--%s->(%s)" % (from_gid, label, to_gid)

    def dump(self):
        return json.dumps({
            "gid": self.gid,
            "label": self.label,
            "from": self.from_gid,
            "to": self.to_gid,
            "data": self.data
        })
