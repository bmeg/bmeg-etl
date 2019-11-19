from bmeg import Allele, Project

genome = 'ncbi_build'
chromosome = 'chromosome'
start = 'start_position'
end = 'end_position'
reference_bases = 'reference_allele'
alternate_bases = 'tumor_seq_allele2'
strand = 'strand'
amino_acids = 'amino_acids'
codons = 'codons'
cds_position = 'cds_position'
cdna_position = 'cdna_position'
protein_position = 'protein_position'
exon_number = 'exon_number'
variant_classification = 'variant_classification'
variant_type = 'variant_type'
biotype = 'biotype'
impact = 'impact'
dbsnp_rs = 'dbsnp_rs'
hgnc_id = 'hgnc_id'
hgvsc = 'hgvsc'
hgvsp = 'hgvsp'
hgvsp_short = 'hgvsp_short'
ensembl_gene = 'gene'
ensembl_transcript = 'transcript_id'
ensembl_protein = 'ensp'
hugo_symbol = 'hugo_symbol'
all_effects = 'all_effects'
polyphen = 'polyphen'
sift = 'sift'
exac_af = 'exac_af'
exac_af_afr = 'exac_af_afr'
exac_af_amr = 'exac_af_amr'
exac_af_adj = 'exac_af_adj'
exac_af_eas = 'exac_af_eas'
exac_af_fin = 'exac_af_fin'
exac_af_nfe = 'exac_af_nfe'
exac_af_oth = 'exac_af_oth'
exac_af_sas = 'exac_af_sas'


def make_minimal_allele(maf_line, alternate_bases="tumor_seq_allele2"):
    assert isinstance(maf_line, dict)
    maf_line = {k.lower(): v for k, v in maf_line.items()}
    alternate_bases = alternate_bases.lower()
    maf_line[genome] = "GRCh37"
    if reference_bases not in maf_line and alternate_bases not in maf_line:
        raise ValueError("ref and alt bases are missing for row: %s", maf_line)
    if reference_bases not in maf_line:
        maf_line[reference_bases] = "."
    if alternate_bases not in maf_line:
        maf_line[alternate_bases] = "."
    return Allele(
        id=Allele.make_gid(
            maf_line[genome], maf_line[chromosome], maf_line[start],
            maf_line[reference_bases], maf_line[alternate_bases],
        ),
        genome=maf_line[genome],
        chromosome=maf_line[chromosome],
        start=int(maf_line[start]),
        reference_bases=maf_line[reference_bases],
        alternate_bases=maf_line[alternate_bases],
        project_id=Project.make_gid('Reference')
    )


def make_allele(maf_line):
    assert isinstance(maf_line, dict)
    maf_line = {k.lower(): v for k, v in maf_line.items()}
    if reference_bases not in maf_line and alternate_bases not in maf_line:
        raise ValueError("ref and alt bases are missing for row: %s", maf_line)
    if reference_bases not in maf_line:
        maf_line[reference_bases] = "."
    if alternate_bases not in maf_line:
        maf_line[alternate_bases] = "."
    return Allele(
        id=Allele.make_gid(
            maf_line[genome], maf_line[chromosome], maf_line[start],
            maf_line[reference_bases], maf_line[alternate_bases],
        ),
        genome=maf_line[genome],
        chromosome=maf_line[chromosome],
        start=int(maf_line[start]),
        end=int(maf_line[end]),
        strand=maf_line.get(strand),
        reference_bases=maf_line[reference_bases],
        alternate_bases=maf_line[alternate_bases],
        amino_acids=maf_line.get(amino_acids),
        codons=maf_line.get(codons),
        cds_position=maf_line.get(cds_position),
        cdna_position=maf_line.get(cdna_position),
        protein_position=maf_line.get(protein_position),
        exon_number=maf_line.get(exon_number),
        variant_classification=maf_line.get(variant_classification),
        variant_type=maf_line.get(variant_type),
        biotype=maf_line.get(biotype),
        impact=maf_line.get(impact),
        dbsnp_rs=maf_line.get(dbsnp_rs),
        hgnc_id=maf_line.get(hgnc_id),
        hgvsc=maf_line.get(hgvsc),
        hgvsp=maf_line.get(hgvsp),
        hgvsp_short=maf_line.get(hgvsp_short),
        ensembl_gene=maf_line.get(ensembl_gene),
        ensembl_transcript=maf_line.get(ensembl_transcript),
        ensembl_protein=maf_line.get(ensembl_protein),
        hugo_symbol=maf_line.get(hugo_symbol),
        all_effects=maf_line.get(all_effects),
        polyphen=maf_line.get(polyphen),
        sift=maf_line.get(sift),
        exac_af=float(maf_line.get(exac_af)) if maf_line.get(exac_af) else None,
        exac_af_adj=float(maf_line.get(exac_af_adj)) if maf_line.get(exac_af_adj) else None,
        exac_af_afr=float(maf_line.get(exac_af_afr)) if maf_line.get(exac_af_afr) else None,
        exac_af_amr=float(maf_line.get(exac_af_amr)) if maf_line.get(exac_af_amr) else None,
        exac_af_eas=float(maf_line.get(exac_af_eas)) if maf_line.get(exac_af_eas) else None,
        exac_af_fin=float(maf_line.get(exac_af_fin)) if maf_line.get(exac_af_fin) else None,
        exac_af_nfe=float(maf_line.get(exac_af_nfe)) if maf_line.get(exac_af_nfe) else None,
        exac_af_oth=float(maf_line.get(exac_af_oth)) if maf_line.get(exac_af_oth) else None,
        exac_af_sas=float(maf_line.get(exac_af_sas)) if maf_line.get(exac_af_sas) else None,
        project_id=Project.make_gid('Reference')
    )


def make_variant_call_data(maf_line, methods):
    assert isinstance(maf_line, dict)
    maf_line = {k.lower(): v for k, v in maf_line.items()}
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
        'reference_allele': 'ref',
        'tumor_seq_allele2': 'alt',
        'filter': 'filter',
    }
    info = {
        "methods": methods
    }
    for k, kn in call_str_keys.items():
        val = maf_line.get(k, None)
        if val == "." or val == "" or val is None:
            val = None
        info[kn] = val
    for k, kn in call_int_keys.items():
        val = maf_line.get(k, None)
        if val == "." or val == "" or val is None:
            info[kn] = None
        else:
            info[kn] = int(val)
    return info
