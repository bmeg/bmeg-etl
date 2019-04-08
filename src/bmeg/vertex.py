from dataclasses import dataclass, asdict
from enum import Enum
import hashlib
from typing import Union
from dacite import from_dict as dacite_from_dict
import re
from bmeg.gid import GID
from bmeg.utils import enforce_types


@enforce_types
@dataclass(frozen=True)
class Vertex:
    def label(self):
        return self.__class__.__name__

    def asdict(self):
        return asdict(self)

    @classmethod
    def from_dict(cls, data):
        if data:
            return dacite_from_dict(data_class=cls, data=data)
        return None


@enforce_types
@dataclass(frozen=True)
class Callset(Vertex):
    tumor_aliquot_id: str
    normal_aliquot_id: Union[None, str]
    source: str

    def gid(self):
        return Callset.make_gid(self.tumor_aliquot_id,
                                self.normal_aliquot_id,
                                self.source)

    @classmethod
    def make_gid(cls, tumor_aliquot_id, normal_aliquot_id, source):
        return GID("%s:%s:%s:%s" % (cls.__name__, source,
                                    tumor_aliquot_id, normal_aliquot_id))


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
    strand: str
    reference_bases: str
    alternate_bases: str
    hugo_symbol: str
    ensembl_transcript: Union[None, str] = None
    type: Union[None, str] = None
    effect: Union[None, str] = None
    dbSNP_RS: Union[None, str] = None

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
class Gene(Vertex):
    gene_id: str
    symbol: str
    description: str
    chromosome: str
    start: int
    end: int
    strand: str
    genome: str

    def gid(self):
        return Gene.make_gid(self.gene_id)

    @classmethod
    def make_gid(cls, gene_id):
        if not gene_id.startswith("ENSG"):
            raise ValueError("not an emsembl gene id: {}".format(gene_id))
        if gene_id.count(".") != 0:
            raise ValueError("version numbers not allowed")
        return GID("%s" % (gene_id))


@enforce_types
@dataclass(frozen=True)
class Transcript(Vertex):
    transcript_id: str
    gene_id: str
    chromosome: str
    start: int
    end: int
    strand: str
    genome: str
    biotype: str

    def gid(self):
        return Transcript.make_gid(self.transcript_id)

    @classmethod
    def make_gid(cls, transcript_id):
        if not transcript_id.startswith("ENST"):
            raise ValueError("not an emsembl transcript id: {}".format(transcript_id))
        if transcript_id.count(".") != 0:
            raise ValueError("version numbers not allowed")
        return GID("%s" % (transcript_id))


@enforce_types
@dataclass(frozen=True)
class Exon(Vertex):
    exon_id: str
    transcript_id: list
    chromosome: str
    start: int
    end: int
    strand: str
    genome: str

    def gid(self):
        return Exon.make_gid(self.exon_id)

    @classmethod
    def make_gid(cls, exon_id):
        if not exon_id.startswith("ENSE"):
            raise ValueError("not an emsembl exon id: {}".format(exon_id))
        if exon_id.count(".") != 0:
            raise ValueError("version numbers not allowed: %s" % (exon_id))
        return GID("%s" % (exon_id))


@enforce_types
@dataclass(frozen=True)
class Protein(Vertex):
    protein_id: str
    transcript_id: str
    genome: str
    uniprot_id: Union[None, str] = None

    def gid(self):
        return Protein.make_gid(self.protein_id)

    @classmethod
    def make_gid(cls, protein_id):
        if not protein_id.startswith("ENSP"):
            raise ValueError("not an emsembl protein id: {}".format(protein_id))
        if protein_id.count(".") != 0:
            raise ValueError("version numbers not allowed")
        return GID("%s" % (protein_id))


@enforce_types
@dataclass(frozen=True)
class ProteinStructure(Vertex):
    protein_id: str

    def gid(self):
        return ProteinStructure.make_gid(self.protein_id)

    @classmethod
    def make_gid(cls, protein_id):
        return GID("%s:%s" % ("PDB", protein_id))


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
        return GID("%s:%s" % (cls.__name__, accession))


class ExpressionMetric(str, Enum):
    TPM = "TPM"
    RPKM = "RKPM"
    GENE_TPM = "GENE_TPM"


@enforce_types
@dataclass(frozen=True)
class Expression(Vertex):
    id: str
    source: str
    metric: ExpressionMetric
    method: str
    values: dict

    def gid(self):
        return Expression.make_gid(self.source, self.id)

    @classmethod
    def make_gid(cls, source, id):
        return GID("%s:%s:%s" % (cls.__name__, source, id))


@enforce_types
@dataclass(frozen=True)
class TranscriptExpression(Expression):

    def gid(self):
        return TranscriptExpression.make_gid(self.source, self.id)

    @classmethod
    def make_gid(cls, source, id):
        return GID("%s:%s:%s" % (cls.__name__, source, id))


@enforce_types
@dataclass(frozen=True)
class GeneExpression(Expression):

    def gid(self):
        return GeneExpression.make_gid(self.source, self.id)

    @classmethod
    def make_gid(cls, source, id):
        return GID("%s:%s:%s" % (cls.__name__, source, id))


@enforce_types
@dataclass(frozen=True)
class PFAMClan(Vertex):
    accession: str
    id: str
    description: str

    def gid(self):
        return PFAMClan.make_gid(self.accession)

    @classmethod
    def make_gid(cls, accession):
        return GID("%s:%s" % (cls.__name__, accession))


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
class Case(Vertex):
    case_id: str
    gdc_attributes: Union[None, dict] = None
    gtex_attributes: Union[None, dict] = None
    ccle_attributes: Union[None, dict] = None

    def gid(self):
        return Case.make_gid(self.case_id)

    @classmethod
    def make_gid(cls, case_id):
        return GID("%s:%s" % (cls.__name__, case_id))


@enforce_types
@dataclass(frozen=True)
class Sample(Vertex):
    sample_id: str
    gdc_attributes: Union[None, dict] = None
    ccle_attributes: Union[None, dict] = None
    gtex_attributes: Union[None, dict] = None

    def gid(self):
        return Sample.make_gid(self.sample_id)

    @classmethod
    def make_gid(cls, sample_id):
        return GID("%s:%s" % (cls.__name__, sample_id))


@enforce_types
@dataclass(frozen=True)
class Aliquot(Vertex):
    aliquot_id: str
    gdc_attributes: Union[None, dict] = None

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
        return GID("%s:%s" % (cls.__name__, go_id))


@enforce_types
@dataclass(frozen=True)
class Program(Vertex):
    program_id: str
    gdc_attributes: Union[None, dict] = None

    def gid(self):
        return Program.make_gid(self.program_id)

    @classmethod
    def make_gid(cls, program_id):
        return GID("%s:%s" % (cls.__name__, program_id))


@enforce_types
@dataclass(frozen=True)
class Project(Vertex):
    project_id: str
    gdc_attributes: Union[None, dict] = None

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
    source: str
    sample_id: str
    compound_id: str

    doses_um: Union[None, list] = None
    activity_data_median: Union[None, list] = None
    activity_sd: Union[None, list] = None
    num_data: Union[None, float] = None
    fit_type: Union[None, str] = None
    ec50: Union[None, float] = None
    ic50: Union[None, float] = None
    amax: Union[None, float] = None
    act_area: Union[None, float] = None

    def gid(self):
        return DrugResponse.make_gid(self.source, self.sample_id, self.compound_id)

    @classmethod
    def make_gid(cls, source, sample_id, compound_id):
        return GID("%s:%s:%s:%s" % (cls.__name__, source, sample_id, compound_id))


@enforce_types
@dataclass(frozen=True)
class Compound(Vertex):
    term_id: str  # compound:CID60823
    term: str  # Atorvastatin
    name: Union[None, str] = None  # some unoffical name e.g. Lipitor, tulip

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
    name: Union[None, str] = None  # some unoffical name e.g. cancer of the breast

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
    url: str  # https://www.ncbi.nlm.nih.gov/pubmed/18451181
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
        rec = re.compile(r"https?://(www\.)?")
        url = rec.sub("", url).strip()
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
class GenomicFeature(Vertex):
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
        return GenomicFeature.make_gid(self.genome, self.chromosome, self.start, self.end, self.annotations, self.myvariantinfo,
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


@enforce_types
@dataclass(frozen=True)
class File(Vertex):
    """ a storage object
    """
    md5: str
    path: str

    def gid(self):
        return File.make_gid(self.md5)

    @classmethod
    def make_gid(cls, md5):
        return GID("%s:%s" % (cls.__name__, md5))


@enforce_types
@dataclass(frozen=True)
class Command(Vertex):
    """ a storage object
    """
    md5: str
    cmd: str
    filename: str

    def gid(self):
        return Command.make_gid(self.md5)

    @classmethod
    def make_gid(cls, md5):
        return GID("%s:%s" % (cls.__name__, md5))


@enforce_types
@dataclass(frozen=True)
class CopyNumberAlteration(Vertex):
    """Gene level copy number alterations
    """
    id: str
    source: str
    method: str
    values: dict

    def gid(self):
        return CopyNumberAlteration.make_gid(self.source, self.id)

    @classmethod
    def make_gid(cls, source, id):
        return GID("%s:%s:%s" % (cls.__name__, source, id))


@enforce_types
@dataclass(frozen=True)
class Methylation(Vertex):
    id: str
    source: str
    metric: str
    method: str
    values: dict

    def gid(self):
        return Methylation.make_gid(self.source, self.id)

    @classmethod
    def make_gid(cls, source, id):
        return GID("%s:%s:%s" % (cls.__name__, source, id))


@enforce_types
@dataclass(frozen=True)
class MethylationProbe(Vertex):
    id: str
    target: Union[None, str]
    chromosome: str
    position: int

    def gid(self):
        return MethylationProbe.make_gid(self.id)

    @classmethod
    def make_gid(cls, id):
        return GID("%s:%s" % (cls.__name__, id))
