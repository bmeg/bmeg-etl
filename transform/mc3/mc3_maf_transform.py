
""" transform a maf file into vertexs[variant, allele]   """
from bmeg.vertex import Callset, Gene
from bmeg.edge import SomaticVariant

from bmeg.maf.maf_transform import get_value, MAFTransformer, maf_default_argument_parser
from bmeg.emitter import new_emitter
from bmeg.util.logging import default_logging
from bmeg.ioutils import reader
import json
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
    # callset source
    SOURCE = 'mc3'
    DEFAULT_PREFIX = SOURCE
    DEFAULT_MAF_FILE = 'source/mc3/mc3.v0.2.8.PUBLIC.maf.gz'

    def barcode_to_aliquot_id(self, barcode):
        """ create tcga sample barcode """
        global ALIQUOT_CONVERSION_TABLE
        if barcode not in ALIQUOT_CONVERSION_TABLE:
            logging.warning('{} not found in ALIQUOT_CONVERSION_TABLE'.format(barcode))
        return ALIQUOT_CONVERSION_TABLE[barcode]

    def create_gene_gid(self, line):  # pragma nocover
        """ override, create gene_gid from line """
        ensembl_id = line.get('Gene', None)
        return Gene.make_gid(gene_id=ensembl_id)

    def allele_call_maker(self, allele, line=None):
        """ create call from line """
        info = {}
        for k, kn in MC3_EXTENSION_CALLSET_KEYS.items():
            info[kn] = get_value(line, k, None)
        for k, kn in MC3_EXTENSION_CALLSET_INT_KEYS.items():
            info[kn] = int(get_value(line, k, None))
        call_methods = line.get('CENTERS', '')
        call_methods = [call_method.replace('*', '') for call_method in call_methods.split("|")]
        info['methods'] = call_methods
        return SomaticVariant(**info)

    def callset_maker(self, allele, source, centerCol, method, line):
        """ create callset from line """
        tumor_aliquot_gid = self.barcode_to_aliquot_id(line['Tumor_Sample_Barcode'])
        normal_aliquot_gid = self.barcode_to_aliquot_id(line['Matched_Norm_Sample_Barcode'])

        sample_callset = Callset(tumor_aliquot_id=tumor_aliquot_gid,
                                 normal_aliquot_id=normal_aliquot_gid,
                                 source=source)
        sample_call = (self.allele_call_maker(allele, line), sample_callset.gid())

        return [sample_call], [sample_callset]


def transform(mafpath, prefix, gdc_aliquot_path, source=MC3_MAFTransformer.SOURCE, emitter_name='json', skip=0, transformer=MC3_MAFTransformer()):
    """ entry point """

    # ensure that we have a lookup from CCLE native barcode to gdc derived uuid
    global ALIQUOT_CONVERSION_TABLE
    with reader(gdc_aliquot_path) as f:
        for line in f:
            aliquot = json.loads(line)
            ALIQUOT_CONVERSION_TABLE[aliquot['data']['gdc_attributes']['submitter_id']] = aliquot['data']['aliquot_id']
    emitter = new_emitter(name=emitter_name, directory=prefix)
    transformer.maf_convert(emitter=emitter, mafpath=mafpath, skip=skip, source=source)

    emitter.close()


def main(transformer=MC3_MAFTransformer()):  # pragma: no cover
    parser = maf_default_argument_parser(transformer)
    parser.add_argument('--gdc_aliquot_path', type=str,
                        help='Path to the directory containing gdc.Aliquot.Vertex.json',
                        default='outputs/gdc/Aliquot.Vertex.json.gz')
    # We don't need the first argument, which is the program name
    options = parser.parse_args(sys.argv[1:])
    default_logging(options.loglevel)
    if not options.maf_file:
        options.maf_file = 'source/mc3/mc3.v0.2.8.PUBLIC.maf.gz'

    transform(mafpath=options.maf_file,
              prefix=options.prefix,
              skip=options.skip,
              source=transformer.SOURCE,
              emitter_name=options.emitter,
              gdc_aliquot_path=options.gdc_aliquot_path,
              transformer=transformer)


if __name__ == '__main__':  # pragma: no cover
    main()
