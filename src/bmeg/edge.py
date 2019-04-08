from dataclasses import dataclass

from bmeg.gid import GID
from bmeg.utils import enforce_types
from dacite import from_dict as dacite_from_dict


@enforce_types
@dataclass(frozen=True)
class Edge:
    def label(self):
        return self.__class__.__name__

    @classmethod
    @enforce_types
    def make_gid(cls, from_gid: GID, to_gid: GID):
        return "(%s)--%s->(%s)" % (from_gid, cls.__name__, to_gid)

    @classmethod
    def from_dict(cls, data):
        if data:
            return dacite_from_dict(data_class=cls, data=data)
        return None


@enforce_types
@dataclass(frozen=True)
class AlleleIn(Edge):
    """
    Allele -> Gene
    """
    def label(self):
        return "alleleIn"


@enforce_types
@dataclass(frozen=True)
class AlleleCall(Edge):
    """
    Callset -> Allele
    """
    ref: str
    alt: str
    t_depth: int
    t_ref_count: int
    t_alt_count: int
    n_depth: int
    n_ref_count: int
    n_alt_count: int
    filter: str
    methods: list
    ensembl_transcript: str
    ensembl_gene: str

    def label(self):
        return "hasSomaticVariant"


@enforce_types
@dataclass(frozen=True)
class MethylationOf(Edge):
    """
    Methylation -> Aliquot
    """

    def label(self):
        return "methylationOf"


@enforce_types
@dataclass(frozen=True)
class MethylationProbeFor(Edge):
    """
    MethylationProbe -> Gene
    """
    def label(self):
        return "methylationProbeOf"


@enforce_types
@dataclass(frozen=True)
class TranscriptFor(Edge):
    """
    Transcript -> Gene
    """
    def label(self):
        return "transcriptOf"


@enforce_types
@dataclass(frozen=True)
class ExonFor(Edge):
    """
    Exon -> Transcript
    """
    def label(self):
        return "exonOf"


@enforce_types
@dataclass(frozen=True)
class ProteinFor(Edge):
    """
    Protein -> Transcript
    """
    def label(self):
        return "proteinOf"


@enforce_types
@dataclass(frozen=True)
class StructureFor(Edge):
    """
    ProteinStructure -> Protein 
    """
    def label(self):
        return "structureOf"


@enforce_types
@dataclass(frozen=True)
class COCAClusterFor(Edge):
    """
    COCACluster -> Case
    """
    def label(self):
        return "hasLabel"


@enforce_types
@dataclass(frozen=True)
class CallsetFor(Edge):
    """
    Callset -> Sample
    """
    def label(self):
        return "callsetOf"    


@enforce_types
@dataclass(frozen=True)
class PFAMClanMember(Edge):
    """
    PFAMClan -> PFAMFamily
    """
    def label(self):
        return "clanOf"


@enforce_types
@dataclass(frozen=True)
class PFAMAlignment(Edge):
    """
    Protein -> PFAMFamily
    """
    start: int
    end: int

    def label(self):
        return "hasFamily"


@enforce_types
@dataclass(frozen=True)
class InProgram(Edge):
    """
    Project -> Program
    """
    def label(self):
        return "projectOf"


@enforce_types
@dataclass(frozen=True)
class InProject(Edge):
    """
    Case -> Project
    """
    def label(self):
        return "caseOf"


@enforce_types
@dataclass(frozen=True)
class SampleFor(Edge):
    """
    Sample -> Case
    """
    def label(self):
        return "sampleOf"


@enforce_types
@dataclass(frozen=True)
class AliquotFor(Edge):
    """
    Aliquot -> Sample
    """
    def label(self):
        return "aliquotOf"


@enforce_types
@dataclass(frozen=True)
class GeneOntologyAnnotation(Edge):
    """
    GenoOntologyTerm -> Gene
    """
    evidence: str
    title: str
    references: list

    def label(self):
        return "annotationOf"


@enforce_types
@dataclass(frozen=True)
class GeneOntologyIsA(Edge):
    """
    GenoOntologyTerm -> GenoOntologyTerm
    """
    def label(self):
        return "subclassOf"


@enforce_types
@dataclass(frozen=True)
class GeneExpressionOf(Edge):
    """
    GeneExpression -> Aliquot
    """
    def label(self):
        return "geneExpressionOf"


@enforce_types
@dataclass(frozen=True)
class TranscriptExpressionOf(Edge):
    """
    TranscriptExpression -> Aliquot
    """
    def label(self):
        return "transcriptExpressionOf"


@enforce_types
@dataclass(frozen=True)
class ResponseIn(Edge):
    """
    DrugResponse -> Aliquot
    """
    def label(self):
        return "responseIn"


@enforce_types
@dataclass(frozen=True)
class ResponseTo(Edge):
    """
    DrugResponse -> Compound
    """
    def label(self):
        return "responseTo"


@enforce_types
@dataclass(frozen=True)
class TestedIn(Edge):
    """
    Compound -> Project
    """
    pass


@enforce_types
@dataclass(frozen=True)
class HasSupportingReference(Edge):
    """
    G2PAssociation -> Publication
    """
    def label(self):
        return "hasSupportingReference"


@enforce_types
@dataclass(frozen=True)
class HasGeneFeature(Edge):
    """
    G2PAssociation -> Gene
    """
    def label(self):
        return "hasGeneFeature"


@enforce_types
@dataclass(frozen=True)
class HasAlleleFeature(Edge):
    """
    G2PAssociation -> Allele
    """
    pass


@enforce_types
@dataclass(frozen=True)
class HasPhenotype(Edge):
    """
    G2PAssociation -> Phenotype
    """
    def label(self):
        return "hasPhenotype"


@enforce_types
@dataclass(frozen=True)
class HasEnvironment(Edge):
    """
    G2PAssociation -> Compound
    """
    def label(self):
        return "hasEnvironment"


@enforce_types
@dataclass(frozen=True)
class HasGenomicFeatureFeature(Edge):
    """
    G2PAssociation -> Gene
    """
    pass


@enforce_types
@dataclass(frozen=True)
class GenomicFeatureIn(Edge):
    """
    GenomicFeature -> Gene
    """
    pass


@enforce_types
@dataclass(frozen=True)
class PhenotypeOf(Edge):
    """
    Aliquot -> Phenotype
    """
    pass


@enforce_types
@dataclass(frozen=True)
class TreatedWith(Edge):
    """
    Case -> Compound
    """
    def label(self):
        return "treatedWith"

@enforce_types
@dataclass(frozen=True)
class Reads(Edge):
    """
    Command -> File
    """
    pass


@enforce_types
@dataclass(frozen=True)
class Writes(Edge):
    """
    Command -> File
    """
    pass


@enforce_types
@dataclass(frozen=True)
class CopyNumberAlterationOf(Edge):
    """
    CopyNumberAlteration -> Aliquot
    """
    def label(self):
        return "copyNumberAlterationOf"
