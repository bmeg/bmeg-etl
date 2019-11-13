import hashlib
import re
from functools import wraps


def default_gid(class_name: str, submitter_id: str):
    return "{}:{}".format(class_name, submitter_id)


def gene_gid(submitter_id: str):
    if not submitter_id.startswith("ENSG"):
        raise ValueError("not an ensembl id: {}".format(submitter_id))
    if submitter_id.count(".") != 0:
        raise ValueError("version numbers not allowed: %s" % (submitter_id))
    return "{}".format(submitter_id)


def transcript_gid(submitter_id: str):
    if not submitter_id.startswith("ENST"):
        raise ValueError("not an ensembl id: {}".format(submitter_id))
    if submitter_id.count(".") != 0:
        raise ValueError("version numbers not allowed: %s" % (submitter_id))
    return "{}".format(submitter_id)


def exon_gid(submitter_id: str):
    if not submitter_id.startswith("ENSE"):
        raise ValueError("not an ensembl id: {}".format(submitter_id))
    if submitter_id.count(".") != 0:
        raise ValueError("version numbers not allowed: %s" % (submitter_id))
    return "{}".format(submitter_id)


def protein_gid(submitter_id: str):
    if not submitter_id.startswith("ENSP"):
        raise ValueError("not an ensembl id: {}".format(submitter_id))
    if submitter_id.count(".") != 0:
        raise ValueError("version numbers not allowed: %s" % (submitter_id))
    return "{}".format(submitter_id)


def drugresponse_gid(source: str, submitter_id: str, compound_id: str):
    return "DrugResponse:{}:{}:{}".format(source, submitter_id, compound_id)


def somatic_callset_gid(source: str, tumor_aliquot_id: str, normal_aliquot_id: str):
    return "SomaticCallset:{}:{}:{}".format(source, tumor_aliquot_id, normal_aliquot_id)


def allele_gid(genome: str, chromosome: str, start: int, end: int,
               reference_bases: str, alternate_bases: str):
    vid = "{}:{}:{}:{}:{}:{}".format(genome, chromosome, start, end, reference_bases, alternate_bases)
    vid = vid.encode('utf-8')
    vidhash = hashlib.sha1()
    vidhash.update(vid)
    vidhash = vidhash.hexdigest()
    return "Allele:{}".format(vidhash)


def genomic_feature_gid(genome: str, chromosome: str, start: int, end: int,
                        type: str, name: str):
    vid = "{}:{}:{}:{}:{}:{}".format(genome, chromosome, start, end, type, name)
    vid = vid.encode('utf-8')
    vidhash = hashlib.sha1()
    vidhash.update(vid)
    vidhash = vidhash.hexdigest()
    return "GenomicFeature:{}".format(vidhash)


def g2p_association_gid(source: str, description: str, evidence_label: str, response_type: str,
                        oncogenic: str, source_document: str, source_url: str):
    a = [p for p in [source, description, evidence_label, response_type, oncogenic, source_document, source_url] if p]
    m = hashlib.sha1()
    m.update(':'.join(a).encode('utf-8'))
    return "G2PAssociation:{}".format(m.hexdigest())


def publication_gid(url: str):
    rec = re.compile(r"https?://(www\.)?")
    url = rec.sub("", url).strip()
    return "Publication:{}".format(url)


def pathway_gid(url: str):
    rec = re.compile(r"https?://(www\.)?")
    url = rec.sub("", url).strip()
    return "{}".format(url)


def pathway_component_gid(url: str):
    rec = re.compile(r"https?://(www\.)?")
    url = rec.sub("", url).strip()
    return "{}".format(url)


def interaction_gid(participant_a: str, interaction_type: str, participant_b: str):
    return "Interaction:{}:{}:{}".format(participant_a, interaction_type, participant_b)


# cast the result of the above gid functions to the proper type
def cast_gid(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        cls = args[0]
        args = args[1:]
        return cls._gid_cls(func(*args, **kwargs))
    return wrapper


# these will be used as the make_gid class method on schema generated vertex classes
gid_factories = {
    'Gene': gene_gid,
    'Transcript': transcript_gid,
    'Protein': protein_gid,
    'Exon': exon_gid,
    'DrugResponse': drugresponse_gid,
    'SomaticCallset': somatic_callset_gid,
    'Allele': allele_gid,
    'GenomicFeature': genomic_feature_gid,
    'G2PAssociation': g2p_association_gid,
    'Publication': publication_gid,
    'Pathway': pathway_gid,
    'PathwayComponent': pathway_component_gid,
    'Interaction': interaction_gid
}
