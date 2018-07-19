import json

from dataclasses import dataclass, field

from bmeg.models.gid_models import GID
from bmeg.models.utils import set_gid


@dataclass(frozen=True)
class Edge:
    gid: GID = field(init=False)
    from_gid: GID
    to_gid: GID

    def __post_init__(self):
        set_gid(self, self.__class__.make_gid(self.label, self.from_gid,
                                              self.to_gid))

    @classmethod
    def make_gid(cls, from_gid, to_gid):
        return "(%s)--%s->(%s)" % (from_gid, cls.__name__, to_gid)

    def dump(self):
        if not self.gid:
            raise ValueError("gid is empty")

        data = dict(self.__dict__)
        del data["gid"]
        del data["from_gid"]
        del data["to"]

        return json.dumps({
            "gid": self.gid,
            "label": self.__class__.__name__,
            "from": self.from_gid,
            "to": self.to_gid,
            "data": data
        })


@dataclass(frozen=True)
class VariantInGene(Edge):
    """
    Variant -> Gene
    """
    pass


@dataclass(frozen=True)
class VariantCall(Edge):
    """
    Variant -> Callset
    """
    info: dict


@dataclass(frozen=True)
class CNASegmentOverlapsGene(Edge):
    """
    CNASegment -> Gene
    """
    pass


@dataclass(frozen=True)
class CNASegmentCall(Edge):
    """
    CNASegment -> Callset
    """
    value: float
    pass


@dataclass(frozen=True)
class CNAValueForGene(Edge):
    """
    Gene -> Callset
    """
    value: float
    pass


@dataclass(frozen=True)
class MethlyationProbeValue(Edge):
    """
    MethylationProbe -> Callset
    """
    value: float
    pass


@dataclass(frozen=True)
class TranscriptFor(Edge):
    """
    Transcript -> Gene
    """
    pass


@dataclass(frozen=True)
class ExonFor(Edge):
    """
    Exon -> Transcript
    """
    pass


@dataclass(frozen=True)
class ProteinFor(Edge):
    """
    Protein -> Transcript
    """
    pass
