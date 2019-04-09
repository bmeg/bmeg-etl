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
class SomaticVariant(Edge):
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
class HasMethylation(Edge):
    """
    Aliquot -> Methylation
    """

    def label(self):
        return "hasMethylation"


@enforce_types
@dataclass(frozen=True)
class HasMethylationProbe(Edge):
    """
    Gene -> MethylationProbe
    """
    def label(self):
        return "hasMethylationProbe"


@enforce_types
@dataclass(frozen=True)
class HasTranscript(Edge):
    """
    Gene -> Transcript
    """
    def label(self):
        return "hasTranscript"


@enforce_types
@dataclass(frozen=True)
class HasExon(Edge):
    """
    Transcript -> Exon
    """
    def label(self):
        return "hasExon"


@enforce_types
@dataclass(frozen=True)
class HasProtein(Edge):
    """
    Transcript -> Protein
    """
    def label(self):
        return "hasProtein"


@enforce_types
@dataclass(frozen=True)
class StructureFor(Edge):
    """
    Protein -> ProteinStructure
    """
    def label(self):
        return "hasStructure"


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
class HasCallset(Edge):
    """
    Sample -> hasCallset
    """
    def label(self):
        return "hasCallset"

    
@enforce_types
@dataclass(frozen=True)
class PFAMClanMember(Edge):
    """
    PFAMFamily -> PFAMClan
    """
    def label(self):
        return "hasClan"


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
class HasProject(Edge):
    """
    Program -> Project
    """
    def label(self):
        return "hasProject"


@enforce_types
@dataclass(frozen=True)
class HasCase(Edge):
    """
    Project -> Case
    """
    def label(self):
        return "hasCase"


@enforce_types
@dataclass(frozen=True)
class HasSample(Edge):
    """
    Case -> Sample
    """
    def label(self):
        return "hasSample"


@enforce_types
@dataclass(frozen=True)
class HasAliquot(Edge):
    """
    Sample -> Aliquot
    """
    def label(self):
        return "hasAliquot"


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
class HasGeneExpression(Edge):
    """
    Aliquot -> GeneExpression
    """
    def label(self):
        return "hasGeneExpression"


@enforce_types
@dataclass(frozen=True)
class HasTranscriptExpression(Edge):
    """
    Aliquot -> TranscriptExpression
    """
    def label(self):
        return "hasTranscriptExpression"


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
    def label(self):
        return "testedIn"

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
class HasGenomicFeature(Edge):
    """
    G2PAssociation -> Gene
    """

    def label(self):
        return "hasGenomicFeature"


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
class HasCopyNumberAlteration(Edge):
    """
    Aliquot -> CopyNumberAlteration
    """
    def label(self):
        return "hasCopyNumberAlteration"
