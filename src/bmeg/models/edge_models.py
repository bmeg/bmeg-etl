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
        set_gid(self, self.__class__.make_gid(self.from_gid, self.to_gid))

    @classmethod
    def make_gid(cls, from_gid, to_gid):
        return "(%s)--%s->(%s)" % (from_gid, cls.__name__, to_gid)

    def dump(self, preserve_nulls=False):
        if not self.gid:
            raise ValueError("gid is empty")

        data = dict(self.__dict__)
        del data["gid"]
        del data["from_gid"]
        del data["to_gid"]

        # delete null values
        if not preserve_nulls:
            remove = [k for k in data if data[k] is None]
            for k in remove:
                del data[k]

        return json.dumps({
            "gid": self.gid,
            "label": self.__class__.__name__,
            "from": self.from_gid,
            "to": self.to_gid,
            "data": data
        })


# TODO deprecate
@dataclass(frozen=True)
class VariantIn(Edge):
    """
    Variant -> Gene
    """
    pass


# TODO deprecate
@dataclass(frozen=True)
class VariantCall(Edge):
    """
    Variant -> Callset
    """
    info: dict


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
