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
class VariantIn(Edge):
    """
    Variant -> Gene
    """
    pass


@enforce_types
@dataclass(frozen=True)
class VariantCall(Edge):
    """
    Variant -> Callset
    """
    info: dict


@enforce_types
@dataclass(frozen=True)
class AlleleIn(Edge):
    """
    Variant -> Gene
    """
    pass


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


@enforce_types
@dataclass(frozen=True)
class MethylationOf(Edge):
    """
    Methylation -> Aliquot
    """
    pass


@enforce_types
@dataclass(frozen=True)
class MethylationProbeFor(Edge):
    """
    MethylationProbe -> Gene
    """
    pass


@enforce_types
@dataclass(frozen=True)
class TranscriptFor(Edge):
    """
    Transcript -> Gene
    """
    pass


@enforce_types
@dataclass(frozen=True)
class ExonFor(Edge):
    """
    Exon -> Transcript
    """
    pass


@enforce_types
@dataclass(frozen=True)
class ProteinFor(Edge):
    """
    Protein -> Transcript
    """
    pass


@enforce_types
@dataclass(frozen=True)
class StructureFor(Edge):
    """
    Protein -> ProteinStructure
    """
    pass


@enforce_types
@dataclass(frozen=True)
class COCAClusterFor(Edge):
    """
    COCACluster -> Individual
    """
    pass


@enforce_types
@dataclass(frozen=True)
class BiosampleFor(Edge):
    """
    Biosample -> Individual
    """
    pass


@enforce_types
@dataclass(frozen=True)
class CallsetFor(Edge):
    """
    Callset -> Biosample
    """
    pass


@enforce_types
@dataclass(frozen=True)
class PFAMClanMember(Edge):
    """
    PFAMClan -> PFAMFamily
    """
    pass


@enforce_types
@dataclass(frozen=True)
class PFAMAlignment(Edge):
    """
    Protein -> PFAMFamily
    """
    start: int
    end: int


class InProject(Edge):
    """
    Individual -> Project
    """
    pass


@enforce_types
@dataclass(frozen=True)
class GeneOntologyAnnotation(Edge):
    """
    GenoOntologyTerm -> Gene
    """
    evidence: str
    title: str
    references: list


@enforce_types
@dataclass(frozen=True)
class GeneOntologyIsA(Edge):
    """
    GenoOntologyTerm -> GenoOntologyTerm
    """
    pass


@enforce_types
@dataclass(frozen=True)
class AliquotFor(Edge):
    """
    Aliquot -> Biosample
    """
    pass


@enforce_types
@dataclass(frozen=True)
class GeneExpressionOf(Edge):
    """
    GeneExpression -> Aliquot
    """
    pass


@enforce_types
@dataclass(frozen=True)
class TranscriptExpressionOf(Edge):
    """
    TranscriptExpression -> Aliquot
    """
    pass


@enforce_types
@dataclass(frozen=True)
class ResponseIn(Edge):
    """
    DrugResponse -> Aliquot
    """
    pass


@enforce_types
@dataclass(frozen=True)
class ResponseTo(Edge):
    """
    DrugResponse -> Compound
    """
    pass


@enforce_types
@dataclass(frozen=True)
class HasSupportingReference(Edge):
    """
    G2PAssociation -> Publication
    """
    pass


@enforce_types
@dataclass(frozen=True)
class HasGeneFeature(Edge):
    """
    G2PAssociation -> Gene
    """
    pass


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
    pass


@enforce_types
@dataclass(frozen=True)
class HasEnvironment(Edge):
    """
    G2PAssociation -> Compound
    """
    pass


@enforce_types
@dataclass(frozen=True)
class HasMinimalAlleleFeature(Edge):
    """
    G2PAssociation -> Gene
    """
    pass


@enforce_types
@dataclass(frozen=True)
class MinimalAlleleIn(Edge):
    """
    MinimalAllele -> Gene
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
    Individual -> Compound
    """
    pass


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
    pass


@enforce_types
@dataclass(frozen=True)
class DerivedFrom(Edge):
    """
    Aliquot -> File
    """
    pass
