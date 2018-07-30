import hashlib

from dataclasses import dataclass, field

from bmeg.gid import GID
from bmeg.utils import set_gid, enforce_types


@enforce_types
@dataclass(frozen=True)
class Vertex:
    gid: GID = field(init=False)


@enforce_types
@dataclass(frozen=True)
class Callset(Vertex):
    tumor_biosample_id: str
    normal_biosample_id: str
    call_method: str
    source: str

    def __post_init__(self):
        set_gid(self, Callset.make_gid(self.source,
                                       self.tumor_biosample_id,
                                       self.normal_biosample_id,
                                       self.call_method))

    @classmethod
    def make_gid(cls, source, tumor_biosample_id, normal_biosample_id,
                 call_method):
        return GID("%s:%s:%s:%s:%s" % (cls.__name__, source,
                                       tumor_biosample_id, normal_biosample_id,
                                       call_method))


@enforce_types
@dataclass(frozen=True)
class Allele(Vertex):
    genome: str
    chromosome: str
    start: int
    end: int
    reference_bases: str
    alternate_bases: str
    myvariantinfo: dict

    def __post_init__(self):
        set_gid(self, Allele.make_gid(self.genome, self.chromosome, self.start,
                                      self.end, self.reference_bases,
                                      self.alternate_bases))

    @classmethod
    def make_gid(cls, genome, chromosome, start, end, reference_bases,
                 alternate_bases):
        # TODO -  figure out better hashing strategy
        vid = "%s:%s:%d:%d:%s:%s" % (genome, chromosome,
                                     start, end, reference_bases,
                                     alternate_bases)
        vidhash = hashlib.sha1()
        vidhash.update(vid)
        vidhash = vidhash.hexdigest()
        return GID("%s:%s" % (cls.__name__, vidhash))


@enforce_types
@dataclass(frozen=True)
class CNASegment(Vertex):
    genome: str
    chromosome: str
    start: int
    end: int

    def __post_init__(self):
        set_gid(self, CNASegment.make_gid(self.genome, self.chromosome,
                                          self.start, self.end))

    @classmethod
    def make_gid(cls, callset_id, genome, chromosome, start, end):
        return GID("%s:%s:%s:%d:%d" % (cls.__name__, callset_id, genome, start,
                                       end))


@enforce_types
@dataclass(frozen=True)
class MethylationProbe(Vertex):
    genome: str
    chromosome: str
    start: int
    end: int
    probe_id: str

    def __post_init__(self):
        set_gid(self, MethylationProbe.make_gid(self.genome, self.chromosome,
                                                self.start, self.end,
                                                self.probe_id))

    @classmethod
    def make_gid(cls, genome, chromosome, start, end, probe_id):
        return GID("%s:%s:%d:%d:%s" % (cls.__name__, genome, start, end,
                                       probe_id))


@enforce_types
@dataclass(frozen=True)
class Gene(Vertex):
    ensembl_id: str
    symbol: str
    description: str
    chromosome: str
    start: int
    end: int
    strand: str
    mygeneinfo: dict

    def __post_init__(self):
        set_gid(self, Gene.make_gid(self.ensembl_id))

    @classmethod
    def make_gid(cls, ensembl_id):
        if not ensembl_id.startswith("ENSG"):
            raise ValueError("not an emsembl gene id")
        return GID("%s:%s" % (cls.__name__, ensembl_id))


@enforce_types
@dataclass(frozen=True)
class Transcript(Vertex):
    ensembl_id: str
    gene_id: str
    chromosome: str
    start: int
    end: int
    strand: str

    def __post_init__(self):
        set_gid(self, Transcript.make_gid(self.ensembl_id))

    @classmethod
    def make_gid(cls, ensembl_id):
        if not ensembl_id.startswith("ENST"):
            raise ValueError("not an emsembl transcript id")
        return GID("%s:%s" % (cls.__name__, ensembl_id))


@enforce_types
@dataclass(frozen=True)
class Exon(Vertex):
    ensembl_id: str
    transcript_id: str
    chromosome: str
    start: int
    end: int
    strand: str

    def __post_init__(self):
        set_gid(self, Exon.make_gid(self.ensembl_id))

    @classmethod
    def make_gid(cls, ensembl_id):
        if not ensembl_id.startswith("ENSE"):
            raise ValueError("not an emsembl exon id")
        return GID("%s:%s" % (cls.__name__, ensembl_id))


@enforce_types
@dataclass(frozen=True)
class Protein(Vertex):
    ensembl_id: str
    transcript_id: str

    def __post_init__(self):
        set_gid(self, Protein.make_gid(self.ensembl_id))

    @classmethod
    def make_gid(cls, ensembl_id):
        if not ensembl_id.startswith("ENSP"):
            raise ValueError("not an emsembl protein id")
        return GID("%s:%s" % (cls.__name__, ensembl_id))


@enforce_types
@dataclass(frozen=True)
class COCACluster(Vertex):
    cluster_id: str

    def __post_init__(self):
        set_gid(self, COCACluster.make_gid(self.cluster_id))

    @classmethod
    def make_gid(cls, cluster_id):
        return GID("%s:%s" % (cls.__name__, cluster_id))


@enforce_types
@dataclass(frozen=True)
class Individual(Vertex):
    individual_id: str
    gdc_attributes: dict

    def __post_init__(self):
        set_gid(self, Individual.make_gid(self.individual_id))

    @classmethod
    def make_gid(cls, individual_id):
        return GID("%s:%s" % (cls.__name__, individual_id))


@enforce_types
@dataclass(frozen=True)
class Biosample(Vertex):
    biosample_id: str
    gdc_attributes: dict

    def __post_init__(self):
        set_gid(self, Individual.make_gid(self.biosample_id))

    @classmethod
    def make_gid(cls, biosample_id):
        return GID("%s:%s" % (cls.__name__, biosample_id))
