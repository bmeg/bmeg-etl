from dataclasses import dataclass, field

from bmeg.gid import GID
from bmeg.utils import set_gid


@dataclass(frozen=True)
class Edge:
    gid: GID = field(init=False)
    from_gid: GID
    to_gid: GID

    def __post_init__(self):
        set_gid(self, self.__class__.make_gid(self.from_gid, self.to_gid))

    @classmethod
    def make_gid(cls, from_gid, to_gid):
        return "(%s)--%s->(%s)" % (from_gid, cls.__name__, to_gid)


@dataclass(frozen=True)
class VariantIn(Edge):
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
class CNASegmentOverlaps(Edge):
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
class GeneCNAValueCall(Edge):
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
class MethlyationProbeFor(Edge):
    """
    MethylationProbe -> Gene
    """
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


@dataclass(frozen=True)
class COCAClusterFor(Edge):
    """
    COCACluster -> Individual
    """
    pass


@dataclass(frozen=True)
class BiosampleFor(Edge):
    """
    Biosample -> Individual
    """
    pass


@dataclass(frozen=True)
class CallsetFor(Edge):
    """
    Callset -> Biosample
    """
    pass
