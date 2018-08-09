from dataclasses import dataclass

from bmeg.gid import GID
from bmeg.utils import model, enforce_types


@model()
class Edge:
    def label(self):
        return self.__class__.__name__

    @classmethod
    @enforce_types
    def make_gid(cls, from_gid: GID, to_gid: GID):
        return "(%s)--%s->(%s)" % (from_gid, cls.__name__, to_gid)


@model()
class VariantIn(Edge):
    """
    Variant -> Gene
    """
    pass


@model()
class VariantCall(Edge):
    """
    Variant -> Callset
    """
    info: dict


@model()
class AlleleIn(Edge):
    """
    Variant -> Gene
    """
    pass


@model()
class AlleleCall(Edge):
    """
    Allele -> Callset
    """
    info: dict = None


@model()
class CNASegmentOverlaps(Edge):
    """
    CNASegment -> Gene
    """
    pass


@model()
class CNASegmentCall(Edge):
    """
    CNASegment -> Callset
    """
    value: float


@model()
class GeneCNAValueCall(Edge):
    """
    Gene -> Callset
    """
    value: float


@model()
class MethlyationProbeValue(Edge):
    """
    MethylationProbe -> Callset
    """
    value: float


@model()
class MethlyationProbeFor(Edge):
    """
    MethylationProbe -> Gene
    """
    pass


@model()
class TranscriptFor(Edge):
    """
    Transcript -> Gene
    """
    pass


@model()
class ExonFor(Edge):
    """
    Exon -> Transcript
    """
    pass


@model()
class ProteinFor(Edge):
    """
    Protein -> Transcript
    """
    pass


@model()
class COCAClusterFor(Edge):
    """
    COCACluster -> Individual
    """
    pass


@model()
class BiosampleFor(Edge):
    """
    Biosample -> Individual
    """
    pass


@model()
class CallsetFor(Edge):
    """
    Callset -> Biosample
    """
    pass


@model()
class InProject(Edge):
    """
    Individual -> Project
    """
    pass


@model()
class AliquotFor(Edge):
    """
    Aliquot -> Biosample
    """
    pass
