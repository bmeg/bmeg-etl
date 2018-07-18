import json

from dataclasses import dataclass, field
from typing import NamedTuple

from gid_models import GID
from utils import set_gid


@dataclass(frozen=True)
class Vertex(object):
    gid: GID = field(init=False)

    def dump(self):
        if not self.gid:
            raise ValueError("gid is empty")

        data = dict(self.__dict__)
        del data["gid"]

        return json.dumps({
            "gid": self.gid,
            "label": self.__class__.__name__,
            "data": data
        })


@dataclass(frozen=True)
class Callset(Vertex):
    tumor_biosample_id: str
    normal_biosample_id: str
    call_method: str

    def __post_init__(self):
        set_gid(self, Callset.make_gid(self.tumor_biosample_id,
                                       self.normal_biosample_id,
                                       self.call_method))

    @classmethod
    def make_gid(cls, tumor_biosample_id, normal_biosample_id, call_method):
        return GID(
            "%s:%s:%s:%s" % (cls.__name__, tumor_biosample_id,
                             normal_biosample_id, call_method)
        )


class Allele(NamedTuple, Vertex):
    genome: str
    chromosome: str
    start: int
    end: int
    reference_bases: str
    alternate_bases: str
    myvariantinfo: dict

    @property
    def gid(self):
        return Allele.make_gid(self.genome, self.chromosome, self.start,
                               self.end, self.reference_bases,
                               self.alternate_bases)

    @classmethod
    def make_gid(cls, genome, chromosome, start, end, reference_bases,
                 alternate_bases):
        # TODO -  figure out hashing strategy
        return GID("%s:%s:%s:%d:%d:%s:%s" % (cls.__name__, genome, chromosome,
                                             start, end, reference_bases,
                                             alternate_bases))


class Gene(NamedTuple, Vertex):
    ensembl_id: str
    symbol: str
    description: str
    chromosome: str
    start: int
    end: int
    strand: str
    mygeneinfo: dict

    @property
    def gid(self):
        return Gene.make_gid(self.ensembl_id)

    @classmethod
    def make_gid(cls, ensembl_id):
        if not ensembl_id.startswith("ENSG"):
            raise ValueError("not an emsembl gene id")
        return GID("%s:%s" % (cls.__name__, ensembl_id))


class Transcript(NamedTuple, Vertex):
    ensembl_id: str
    gene_id: str
    chromosome: str
    start: int
    end: int
    strand: str

    @property
    def gid(self):
        return Transcript.make_gid(self.ensembl_id)

    @classmethod
    def make_gid(cls, ensembl_id):
        if not ensembl_id.startswith("ENST"):
            raise ValueError("not an emsembl transcript id")
        return GID("%s:%s" % (cls.__name__, ensembl_id))


class Exon(NamedTuple, Vertex):
    ensembl_id: str
    transcript_id: str
    chromosome: str
    start: int
    end: int
    strand: str

    @property
    def gid(self):
        return Exon.make_gid(self.ensembl_id)

    @classmethod
    def make_gid(cls, ensembl_id):
        if not ensembl_id.startswith("ENSE"):
            raise ValueError("not an emsembl exon id")
        return GID("%s:%s" % (cls.__name__, ensembl_id))


class Protein(NamedTuple, Vertex):
    ensembl_id: str
    transcript_id: str

    @property
    def gid(self):
        return Protein.make_gid(self.ensembl_id)

    @classmethod
    def make_gid(cls, ensembl_id):
        if not ensembl_id.startswith("ENSP"):
            raise ValueError("not an emsembl protein id")
        return GID("%s:%s" % (cls.__name__, ensembl_id))
