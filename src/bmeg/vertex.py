from dataclasses import dataclass, field
from enum import Enum
import hashlib
from typing import Union

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
    normal_biosample_id: Union[None, str]
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
@dataclass(frozen=False)  # note: mutable class
class AlleleAnnotations:
    maf: Union[None, dict] = None   # annotations from standard maf https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/
    mc3: Union[None, dict] = None   # annotations from mc3 maf extension
    ccle: Union[None, dict] = None  # annotations from ccle maf  extension
    myvariantinfo: Union[None, dict] = None  # annotations from biothings


@enforce_types
@dataclass(frozen=True)
class Allele(Vertex):
    genome: str
    chromosome: str
    start: int
    end: int
    reference_bases: str
    alternate_bases: str
    annotations: AlleleAnnotations
    type: Union[None, str] = None
    effect: Union[None, str] = None

    def gid(self):
        return Allele.make_gid(self.genome, self.chromosome, self.start,
                               self.end, self.reference_bases,
                               self.alternate_bases)

    @classmethod
    def from_dict(cls, data):
        data['annotations'] = AlleleAnnotations(**data['annotations'])
        return Allele(**data)

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
class PFAMClan(Vertex):
    accession: str

    def gid(self):
        return PFAMClan.make_gid(self.accession)

    @classmethod
    def make_gid(cls, accession):
        return GID("%s:%s" % ("PFAMCLAN", accession))


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
    gtex_attributes: dict

    def gid(self):
        return Individual.make_gid(self.individual_id)

    @classmethod
    def make_gid(cls, individual_id):
        return GID("%s:%s" % (cls.__name__, individual_id))


@enforce_types
@dataclass(frozen=True)
class Biosample(Vertex):
    biosample_id: str
    # TODO normalize style w/ Allele annotations
    gdc_attributes: dict = field(default_factory=dict)
    ccle_attributes: dict = field(default_factory=dict)
    gtex_attributes: dict = field(default_factory=dict)

    def gid(self):
        return Biosample.make_gid(self.biosample_id)

    @classmethod
    def make_gid(cls, biosample_id):
        return GID("%s:%s" % (cls.__name__, biosample_id))


@enforce_types
@dataclass(frozen=True)
class Aliquot(Vertex):
    aliquot_id: str

    def gid(self):
        return Aliquot.make_gid(self.aliquot_id)

    @classmethod
    def make_gid(cls, aliquot_id):
        return GID("%s:%s" % (cls.__name__, aliquot_id))


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


@enforce_types
@dataclass(frozen=True)
class Project(Vertex):
    project_id: str
    gdc_attributes: dict

    def gid(self):
        return Project.make_gid(self.project_id)

    @classmethod
    def make_gid(cls, project_id):
        return GID("%s:%s" % (cls.__name__, project_id))


class DrugResponseMetric(str, Enum):
    AUC = "AUC"
    IC50 = "IC50"


@enforce_types
@dataclass(frozen=True)
class DrugResponse(Vertex):
    compound_name: str
    sample_id: str
    metric: DrugResponseMetric
    value: Union[None, float]

    def gid(self):
        return DrugResponse.make_gid(self.compound_name, self.sample_id)

    @classmethod
    def make_gid(cls, compound_name, sample_id):
        return GID("%s:%s:%s" % (cls.__name__, compound_name, sample_id))


@enforce_types
@dataclass(frozen=True)
class Compound(Vertex):
    term_id: str  # compound:CID60823
    term: str  # Atorvastatin

    def gid(self):
        return Compound.make_gid(self.term_id)

    @classmethod
    def make_gid(cls, term_id):
        return GID("%s:%s" % (cls.__name__, term_id))


@enforce_types
@dataclass(frozen=True)
class Phenotype(Vertex):
    term_id: str  # DOID:11198 ; MONDO:0007254
    term: str  # DiGeorge syndrome ; breast cancer

    def gid(self):
        return Phenotype.make_gid(self.term_id)

    @classmethod
    def make_gid(cls, term_id):
        return GID("%s:%s" % (cls.__name__, term_id))


@enforce_types
@dataclass(frozen=True)
class G2PAssociation(Vertex):
    source: str  # civic, cgi
    description: str  # asprin cures headaches
    evidence_label: Union[None, str]  # evidence
    response_type: Union[None, str]   # evidence
    oncogenic: Union[None, str]  # for non drug evidence source:[oncokb, brca]
    source_document: Union[None, str]  # stringified document from source
    source_url: Union[None, str]  # link back to original document

    def gid(self):
        return G2PAssociation.make_gid(self.source,
                                       self.description,
                                       self.evidence_label,
                                       self.response_type,
                                       self.oncogenic,
                                       self.source_document,
                                       self.source_url)

    @classmethod
    def make_gid(cls, source, description, evidence_label, response_type,
                 oncogenic, source_document, source_url):
        a = [p for p in [source, description, evidence_label, response_type, oncogenic, source_document, source_url] if p]
        m = hashlib.sha1()
        m.update(':'.join(a).encode('utf-8'))
        return GID("%s:%s" % (cls.__name__, m.hexdigest()))


@enforce_types
@dataclass(frozen=True)
class Publication(Vertex):
    url: str  # http://www.ncbi.nlm.nih.gov/pubmed/18451181
    title: Union[None, str]
    abstract: Union[None, str]
    text: Union[None, str]
    date: Union[None, str]

    author: Union[None, list]
    citation: Union[None, list]

    def gid(self):
        return Publication.make_gid(self.url)

    @classmethod
    def make_gid(cls, url):
        return GID("%s:%s" % (cls.__name__, url))


@enforce_types
@dataclass(frozen=True)
class Deadletter(Vertex):
    """ standard way to log missing data """
    target_label: str  # desired vertex label
    data: dict  # data from source

    def gid(self):
        return Deadletter.make_gid(self.target_label, self.data)

    @classmethod
    def make_gid(cls, target_label, data):
        # create hash from data dict
        data = '%s:%s' % (target_label, str(data))
        datahash = hashlib.sha1()
        datahash.update(data.encode('utf-8'))
        datahash = datahash.hexdigest()
        return GID("%s:%s" % (cls.__name__, datahash))


@enforce_types
@dataclass(frozen=True)
class MinimalAllele(Vertex):
    """ consensus set of minimal variant level data (MVLD)
        inspired by https://www.ncbi.nlm.nih.gov/pubmed/27814769
    """
    genome: Union[None, str] = None
    chromosome: Union[None, str] = None
    start: Union[None, int] = None
    end: Union[None, int] = None
    annotations: Union[None, list] = None
    myvariantinfo: Union[None, dict] = None
    type: Union[None, str] = None
    effect: Union[None, str] = None
    name: str = None

    def gid(self):
        return MinimalAllele.make_gid(self.genome, self.chromosome, self.start, self.end, self.annotations, self.myvariantinfo,
                                      self.type, self.effect, self.name)

    @classmethod
    def make_gid(cls, genome, chromosome, start, end, annotations, myvariantinfo, type, effect, name):
        # TODO
        # figure out better hashing strategy
        vid = "%s:%s:%s:%s:%s:%s:%s:%s:%s" % (genome, chromosome, start, end, annotations, myvariantinfo, type, effect, name)
        vid = vid.encode('utf-8')
        vidhash = hashlib.sha1()
        vidhash.update(vid)
        vidhash = vidhash.hexdigest()
        return GID("%s:%s" % (cls.__name__, vidhash))
