from bmeg import Allele, Project

genome = 'NCBI_Build'
chromosome = 'Chromosome'
start = 'Start_Position'
end = 'End_Position'
reference_bases = 'Reference_Allele'
alternate_bases = 'Tumor_Seq_Allele2'
strand = 'Strand'
amino_acids = 'Amino_acids'
codons = 'Codons'
cds_position = 'CDS_Position'
cdna_position = 'cDNA_Position'
protein_position = 'Protein_Position'
exon_number = 'Exon_Number'
variant_classification = 'Variant_Classification'
variant_type = 'Variant_Type'
biotype = 'BIOTYPE'
impact = 'IMPACT'
dbsnp_rs = 'dbSNP_RS'
hgnc_id = 'HGNC_ID'
hgvsc = 'HGVSc'
hgvsp = 'HGVSp'
hgvsp_short = 'HGVSp_Short'
ensembl_gene = 'Gene'
ensembl_transcript = 'Transcript_ID'
ensembl_protein = 'ENSP'
hugo_symbol = 'Hugo_Symbol'
all_effects = 'all_effects'
polyphen = 'PolyPhen'
sift = 'SIFT'
exac_af = 'ExAC_AF'
exac_af_afr = 'ExAC_AF_AFR'
exac_af_amr = 'ExAC_AF_AMR'
exac_af_adj = 'ExAC_AF_ADJ'
exac_af_eas = 'ExAC_AF_EAS'
exac_af_fin = 'ExAC_AF_FIN'
exac_af_nfe = 'ExAC_AF_NFE'
exac_af_oth = 'ExAC_AF_OTH'
exac_af_sas = 'ExAC_AF_SAS'


def make_minimal_allele(maf_line):
    return Allele(
        id=Allele.make_gid(
            maf_line.get(genome), maf_line.get(chromosome), maf_line.get(start),
            maf_line.get(reference_bases), maf_line.get(alternate_bases),
        ),
        chromosome=maf_line.get(chromosome),
        start=int(maf_line.get(start)),
        reference_bases=maf_line.get(reference_bases),
        alternate_bases=maf_line.get(alternate_bases),
        project_id=Project.make_gid('Reference')
    )


def make_allele(maf_line):
    return Allele(
        id=Allele.make_gid(
            maf_line.get(genome), maf_line.get(chromosome), maf_line.get(start),
            maf_line.get(reference_bases), maf_line.get(alternate_bases),
        ),
        genome=maf_line.get(genome),
        chromosome=maf_line.get(chromosome),
        start=int(maf_line.get(start)),
        end=int(maf_line.get(end)),
        strand=maf_line.get(strand),
        reference_bases=maf_line.get(reference_bases),
        alternate_bases=maf_line.get(alternate_bases),
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
