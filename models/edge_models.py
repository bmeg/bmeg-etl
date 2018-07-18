import json

from gid_models import GID


class GenericEdge(object):
    def __init__(self, from_gid, to_gid):
        if not isinstance(from_gid, GID):
            raise TypeError(
                "expected GID not %s" % (from_gid.__class__.__name__)
            )
        if not isinstance(to_gid, GID):
            raise TypeError(
                "expected GID not %s" % (to_gid.__class__.__name__)
            )

        self.label = self.__class__.__name__
        self.id = "(%s:%s)--%s->(%s:%s)" % (from_gid, self.label, to_gid)
        self.from_gid = from_gid
        self.to_gid = to_gid

    def dump(self):
        return json.dumps({
            "gid": self.id,
            "label": self.label,
            "from": self.from_gid,
            "to": self.to_gid,
            "data": self.data
        })
