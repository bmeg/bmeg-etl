import hashlib

from dataclasses import dataclass

from bmeg.gid import GID
from bmeg.utils import enforce_types


@enforce_types
@dataclass(frozen=True)
class Vertex:
    def label(self):
        return self.__class__.__name__


@enforce_types
@dataclass(frozen=True)
class Callset(Vertex):
    tumor_biosample_id: str
    normal_biosample_id: str
    call_method: str
    source: str

    def gid(self):
        return Callset.make_gid(self.tumor_biosample_id,
                                self.normal_biosample_id,
                                self.call_method,
                                self.source)

    @classmethod
    def make_gid(cls, tumor_biosample_id, normal_biosample_id,
                 call_method, source):
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
    annotations: list = None
    myvariantinfo: dict = None

    def gid(self):
        return Allele.make_gid(self.genome, self.chromosome, self.start,
                               self.end, self.reference_bases,
                               self.alternate_bases)

    @classmethod
    def make_gid(cls, genome, chromosome, start, end, reference_bases,
                 alternate_bases):
        # TODO
        # figure out better hashing strategy
        vid = "%s:%s:%d:%d:%s:%s" % (genome, chromosome,
                                     start, end, reference_bases,
                                     alternate_bases)
        vid = vid.encode('utf-8')
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

    def gid(self):
        return CNASegment.make_gid(self.genome, self.chromosome,
                                   self.start, self.end)

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

    def gid(self):
        return MethylationProbe.make_gid(self.genome, self.chromosome,
                                         self.start, self.end,
                                         self.probe_id)

    @classmethod
    def make_gid(cls, genome, chromosome, start, end, probe_id):
        return GID("%s:%s:%d:%d:%s" % (cls.__name__, genome, start, end,
                                       probe_id))


@enforce_types
@dataclass(frozen=True)
class Gene(Vertex):
    gene_id: str
    symbol: str
    description: str
    chromosome: str
    start: int
    end: int
    strand: str
    mygeneinfo: dict

    def gid(self):
        return Gene.make_gid(self.gene_id)

    @classmethod
    def make_gid(cls, gene_id):
        if not gene_id.startswith("ENSG"):
            raise ValueError("not an emsembl gene id")
        return GID("%s:%s" % (cls.__name__, gene_id))


@enforce_types
@dataclass(frozen=True)
class Transcript(Vertex):
    transcript_id: str
    gene_id: str
    chromosome: str
    start: int
    end: int
    strand: str

    def gid(self):
        return Transcript.make_gid(self.transcript_id)

    @classmethod
    def make_gid(cls, transcript_id):
        if not transcript_id.startswith("ENST"):
            raise ValueError("not an emsembl transcript id")
        return GID("%s:%s" % (cls.__name__, transcript_id))


@enforce_types
@dataclass(frozen=True)
class Exon(Vertex):
    exon_id: str
    transcript_id: str
    chromosome: str
    start: int
    end: int
    strand: str

    def gid(self):
        return Exon.make_gid(self.exon_id)

    @classmethod
    def make_gid(cls, exon_id):
        if not exon_id.startswith("ENSE"):
            raise ValueError("not an emsembl exon id")
        return GID("%s:%s" % (cls.__name__, exon_id))


@enforce_types
@dataclass(frozen=True)
class Protein(Vertex):
    protein_id: str
    transcript_id: str

    def gid(self):
        return Protein.make_gid(self.protein_id)

    @classmethod
    def make_gid(cls, protein_id):
        if not protein_id.startswith("ENSP"):
            raise ValueError("not an emsembl protein id")
        return GID("%s:%s" % (cls.__name__, protein_id))

@enforce_types
@dataclass(frozen=True)
class PFAMFamily(Vertex):
    pfam_id: str
    accession: str
    type: str
    description: str
    comments: str

    def gid(self):
        return PFAMFamily.make_gid(self.accession)

    @classmethod
    def make_gid(cls, accession):
        return GID("%s:%s" % ("PFAM", accession))


@enforce_types
@dataclass(frozen=True)
class COCACluster(Vertex):
    cluster_id: str

    def gid(self):
        return COCACluster.make_gid(self.cluster_id)

    @classmethod
    def make_gid(cls, cluster_id):
        return GID("%s:%s" % (cls.__name__, cluster_id))


@enforce_types
@dataclass(frozen=True)
class Individual(Vertex):
    individual_id: str
    gdc_attributes: dict

    def gid(self):
        return Individual.make_gid(self.individual_id)

    @classmethod
    def make_gid(cls, individual_id):
        return GID("%s:%s" % (cls.__name__, individual_id))


@enforce_types
@dataclass(frozen=True)
class Biosample(Vertex):
    biosample_id: str
    gdc_attributes: dict

    def gid(self):
        return Individual.make_gid(self.biosample_id)

    @classmethod
    def make_gid(cls, biosample_id):
        return GID("%s:%s" % (cls.__name__, biosample_id))


@enforce_types
@dataclass(frozen=True)
class GeneOntologyTerm(Vertex):
    go_id: str
    name: str
    namespace: str
    definition: str
    synonym: list
    xref: list

    def gid(self):
        return GeneOntologyTerm.make_gid(self.go_id)

    @classmethod
    def make_gid(cls, go_id):
        if go_id.startswith("GO:"):
            return GID(go_id)
        return GID("GO:%s" % (go_id))
