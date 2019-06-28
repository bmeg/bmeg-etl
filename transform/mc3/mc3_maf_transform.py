
""" transform a maf file into vertexs[variant, allele]   """
from bmeg import Callset, Gene, Project

from bmeg.maf.maf_transform import get_value, MAFTransformer, maf_default_argument_parser
from bmeg.emitter import new_emitter
from bmeg.util.logging import default_logging
import bmeg.ioutils
import sys
import logging

TUMOR_SAMPLE_BARCODE = "Tumor_Sample_Barcode"  # 15
NORMAL_SAMPLE_BARCODE = "Matched_Norm_Sample_Barcode"  # 16

MC3_EXTENSION_MAF_KEYS = [
    'Verification_Status', 'Validation_Status', 'Mutation_Status', 'Sequencing_Phase', 'Sequence_Source', 'Validation_Method', 'Score', 'BAM_File', 'Sequencer', 'Tumor_Sample_UUID', 'Matched_Norm_Sample_UUID', 'HGVSc', 'HGVSp', 'HGVSp_Short', 'Transcript_ID', 'Exon_Number', 'Allele', 'Gene', 'Feature', 'Feature_type', 'Consequence', 'cDNA_position', 'CDS_position', 'Protein_position', 'Amino_acids', 'Codons', 'Existing_variation', 'ALLELE_NUM', 'DISTANCE', 'STRAND',
    'SYMBOL', 'SYMBOL_SOURCE', 'HGNC_ID', 'BIOTYPE', 'CANONICAL', 'CCDS', 'ENSP', 'SWISSPROT', 'TREMBL', 'UNIPARC', 'RefSeq', 'SIFT', 'PolyPhen', 'EXON', 'INTRON', 'DOMAINS', 'GMAF', 'AFR_MAF', 'AMR_MAF', 'ASN_MAF', 'EAS_MAF', 'EUR_MAF', 'SAS_MAF', 'AA_MAF', 'EA_MAF', 'CLIN_SIG', 'SOMATIC', 'PUBMED', 'MOTIF_NAME', 'MOTIF_POS', 'HIGH_INF_POS', 'MOTIF_SCORE_CHANGE', 'IMPACT', 'PICK', 'VARIANT_CLASS', 'TSL', 'HGVS_OFFSET', 'PHENO', 'MINIMISED', 'ExAC_AF', 'ExAC_AF_AFR', 'ExAC_AF_AMR', 'ExAC_AF_EAS', 'ExAC_AF_FIN', 'ExAC_AF_NFE', 'ExAC_AF_OTH', 'ExAC_AF_SAS', 'GENE_PHENO', 'FILTER', 'COSMIC', 'CENTERS', 'CONTEXT', 'DBVS', 'NCALLERS'
]

MC3_EXTENSION_CALLSET_INT_KEYS = {
    't_depth': 't_depth',
    't_ref_count': 't_ref_count',
    't_alt_count': 't_alt_count',
    'n_depth': 'n_depth',
    'n_ref_count': 'n_ref_count',
    'n_alt_count': 'n_alt_count'
}

MC3_EXTENSION_CALLSET_KEYS = {
    'Reference_Allele': 'ref',
    'Tumor_Seq_Allele1': 'alt',
    'FILTER': 'filter',
    'Gene': 'ensembl_gene',
    'Transcript_ID': 'ensembl_transcript'
}

ALIQUOT_CONVERSION_TABLE = {}


class MC3_MAFTransformer(MAFTransformer):

    # override the alt allele column
    TUMOR_ALLELE = "Tumor_Seq_Allele2"  # 11
    DEFAULT_MAF_FILE = 'source/mc3/mc3.v0.2.8.PUBLIC.maf.gz'

    def barcode_to_aliquot_id(self, barcode):
        """ create tcga sample barcode """
        global ID_CONVERSION_TABLE
        if barcode not in ID_CONVERSION_TABLE:
            logging.warning('{} not found in ID_CONVERSION_TABLE'.format(barcode))
            return barcode
        return ID_CONVERSION_TABLE[barcode]

    def create_gene_gid(self, line):  # pragma nocover
        """ override, create gene_gid from line """
        ensembl_id = line.get('Gene', None)
        return Gene.make_gid(ensembl_id)

    def allele_call_maker(self, line, method):
        """ create call from line """
        info = {}
        for k, kn in MC3_EXTENSION_CALLSET_KEYS.items():
            info[kn] = get_value(line, k, None)
        for k, kn in MC3_EXTENSION_CALLSET_INT_KEYS.items():
            info[kn] = int(get_value(line, k, None))
        call_methods = line.get('CENTERS', '')
        call_methods = [call_method.replace('*', '') for call_method in call_methods.split("|")]
        info['methods'] = call_methods
        return info

    def callset_maker(self, line, method):
        """ create callset from line """
        global PROJECT_CONVERSION_TABLE
        tumor_aliquot_gid = self.barcode_to_aliquot_id(line['Tumor_Sample_Barcode'])
        normal_aliquot_gid = self.barcode_to_aliquot_id(line['Matched_Norm_Sample_Barcode'])
        project_id = PROJECT_CONVERSION_TABLE.get(tumor_aliquot_gid, None)
        sample_callset = Callset(
            id=Callset.make_gid("MC3", tumor_aliquot_gid, normal_aliquot_gid),
            tumor_aliquot_id=tumor_aliquot_gid,
            normal_aliquot_id=normal_aliquot_gid,
            project_id=Project.make_gid(project_id)
        )
        sample_call = (self.allele_call_maker(line, method), sample_callset.gid())
        return [sample_call], [sample_callset]


def transform(mafpath='source/mc3/mc3.v0.2.8.PUBLIC.maf.gz',
              id_lookup_path='source/gdc/id_lookup.tsv',
              project_lookup_path='source/gdc/project_lookup.tsv',
              skip=0,
              emitter_name='json',
              emitter_directory='mc3'):

    # ensure that we have a lookup from TCGA barcode to gdc uuid
    global ID_CONVERSION_TABLE
    ID_CONVERSION_TABLE = bmeg.ioutils.read_lookup(id_lookup_path)

    # ensure that we have a lookup from TCGA id to project id
    global PROJECT_CONVERSION_TABLE
    PROJECT_CONVERSION_TABLE = bmeg.ioutils.read_lookup(project_lookup_path)

    transformer = MC3_MAFTransformer()
    emitter = new_emitter(name=emitter_name, directory=emitter_directory)

    transformer.maf_convert(
        emitter=emitter,
        mafpath=mafpath,
        skip=skip
    )

    emitter.close()


def main(transformer=MC3_MAFTransformer()):  # pragma: no cover
    parser = maf_default_argument_parser(transformer)
    # We don't need the first argument, which is the program name
    options = parser.parse_args(sys.argv[1:])
    default_logging(options.loglevel)
    if not options.maf_file:
        options.maf_file = 'source/mc3/mc3.v0.2.8.PUBLIC.maf.gz'

    transform(
        mafpath=options.maf_file,
        emitter_name=options.emitter,
        skip=options.skip
    )


if __name__ == '__main__':  # pragma: no cover
    main()
