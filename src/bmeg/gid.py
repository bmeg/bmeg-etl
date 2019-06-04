from functools import partial, wraps


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
    'Program': partial(default_gid, 'Program'),
    'Project': partial(default_gid, 'Project'),
    'Case': partial(default_gid, 'Case'),
    'Sample': partial(default_gid, 'Sample'),
    'Aliquot': partial(default_gid, 'Aliquot'),
    'Gene': gene_gid,
    'Transcript': transcript_gid,
    'Protein': protein_gid,
    'Exon': exon_gid,
}
