import hashlib

from bmeg.gid import GID
from bmeg.utils import model, Vertex


@model("tumor_biosample_id normal_biosample_id call_method source")
class Callset(Vertex):
    tumor_biosample_id: str
    normal_biosample_id: str
    call_method: str
    source: str


@model
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


@model("genome chromosome start end")
class CNASegment(Vertex):
    genome: str
    chromosome: str
    start: int
    end: int


@model("genome chromosome start end probe_id")
class MethylationProbe(Vertex):
    genome: str
    chromosome: str
    start: int
    end: int
    probe_id: str


@model("gene_id")
class Gene(Vertex):
    gene_id: str
    symbol: str
    description: str
    chromosome: str
    start: int
    end: int
    strand: str
    mygeneinfo: dict


@model("transcript_id")
class Transcript(Vertex):
    transcript_id: str
    gene_id: str
    chromosome: str
    start: int
    end: int
    strand: str


@model("exon_id")
class Exon(Vertex):
    exon_id: str
    transcript_id: str
    chromosome: str
    start: int
    end: int
    strand: str


@model("protein_id")
class Protein(Vertex):
    protein_id: str
    transcript_id: str


@model("cluster_id")
class COCACluster(Vertex):
    cluster_id: str


@model("individual_id")
class Individual(Vertex):
    individual_id: str
    gdc_attributes: dict


@model("biosample_id")
class Biosample(Vertex):
    biosample_id: str
    gdc_attributes: dict


@model("project_id")
class Project(Vertex):
    project_id: str
    gdc_attributes: dict


@model("aliquot_id")
class Aliquot(Vertex):
    aliquot_id: str
    gdc_attributes: dict
