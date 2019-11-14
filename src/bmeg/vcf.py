from bmeg import Allele, Project

genome = 'NCBI_Build'
chromosome = 'Chromosome'
start = 'Start_Position'
end = 'End_Position'
reference_bases = 'Reference_Allele'
alternate_bases = 'Tumor_Seq_Allele2'
strand = 'Strand'


def make_minimal_allele(vcf_line):
    assert isinstance(vcf_line, dict)
    return Allele(
        id=Allele.make_gid(
            vcf_line.get(genome), vcf_line.get(chromosome), vcf_line.get(start),
            vcf_line.get(reference_bases), vcf_line.get(alternate_bases),
        ),
        chromosome=vcf_line.get(chromosome),
        start=int(vcf_line.get(start)),
        reference_bases=vcf_line.get(reference_bases),
        alternate_bases=vcf_line.get(alternate_bases),
        project_id=Project.make_gid('Reference')
    )


def make_variant_call_data(vcf_line, methods):
    assert isinstance(vcf_line, dict)
    if isinstance(methods, str):
        methods = [methods]
    elif isinstance(methods, list):
        pass
    else:
        raise TypeError("expected a str or list")
    call_int_keys = {
        't_depth': 't_depth',
        't_ref_count': 't_ref_count',
        't_alt_count': 't_alt_count',
        'n_depth': 'n_depth',
        'n_ref_count': 'n_ref_count',
        'n_alt_count': 'n_alt_count'
    }
    call_str_keys = {
        'Reference_Allele': 'ref',
        'Tumor_Seq_Allele2': 'alt',
        'FILTER': 'filter',
    }
    info = {
        "methods": methods
    }
    for k, kn in call_str_keys:
        val = vcf_line.get(k, None)
        if val == "." or val == "" or val is None:
            val = None
        info[kn] = val
    for k, kn in call_int_keys.items():
        val = vcf_line.get(k, None)
        if val == "." or val == "" or val is None:
            info[kn] = None
        else:
            info[kn] = int(val)
    return info
