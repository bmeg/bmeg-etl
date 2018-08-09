from dataclasses import dataclass

from bmeg.gid import GID
from bmeg.utils import enforce_types


@enforce_types
@dataclass(frozen=True)
class Edge:
    def label(self):
        return self.__class__.__name__

    @classmethod
    @enforce_types
    def make_gid(cls, from_gid: GID, to_gid: GID):
        return "(%s)--%s->(%s)" % (from_gid, cls.__name__, to_gid)


@enforce_types
@dataclass(frozen=True)
class VariantIn(Edge):
    """
    Variant -> Gene
    """
    pass


@enforce_types
@dataclass(frozen=True)
class VariantCall(Edge):
    """
    Variant -> Callset
    """
    info: dict


@enforce_types
@dataclass(frozen=True)
class AlleleIn(Edge):
    """
    Variant -> Gene
    """
    pass


@dataclass(frozen=True)
class AlleleCall(Edge):
    """
    Allele -> Callset
    """
    info: dict = None


@dataclass(frozen=True)
class CNASegmentOverlaps(Edge):
    """
    CNASegment -> Gene
    """
    pass


@enforce_types
@dataclass(frozen=True)
class CNASegmentCall(Edge):
    """
    CNASegment -> Callset
    """
    value: float


@enforce_types
@dataclass(frozen=True)
class GeneCNAValueCall(Edge):
    """
    Gene -> Callset
    """
    value: float


@enforce_types
@dataclass(frozen=True)
class MethlyationProbeValue(Edge):
    """
    MethylationProbe -> Callset
    """
    value: float


@enforce_types
@dataclass(frozen=True)
class MethlyationProbeFor(Edge):
    """
    MethylationProbe -> Gene
    """
    pass


@enforce_types
@dataclass(frozen=True)
class TranscriptFor(Edge):
    """
    Transcript -> Gene
    """
    pass


@enforce_types
@dataclass(frozen=True)
class ExonFor(Edge):
    """
    Exon -> Transcript
    """
    pass


@enforce_types
@dataclass(frozen=True)
class ProteinFor(Edge):
    """
    Protein -> Transcript
    """
    pass


@enforce_types
@dataclass(frozen=True)
class COCAClusterFor(Edge):
    """
    COCACluster -> Individual
    """
    pass


@enforce_types
@dataclass(frozen=True)
class BiosampleFor(Edge):
    """
    Biosample -> Individual
    """
    pass


@enforce_types
@dataclass(frozen=True)
class CallsetFor(Edge):
    """
    Callset -> Biosample
    """
    pass


@enforce_types
@dataclass(frozen=True)
class InProject(Edge):
    """
    Individual -> Project
    """
    pass


@enforce_types
@dataclass(frozen=True)
class BiosampleFor(Edge):
    """
    Biosample -> Individual
    """
    pass


@enforce_types
@dataclass(frozen=True)
class AliquotFor(Edge):
    """
    Aliquot -> Biosample
    """
    pass
