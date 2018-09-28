
""" transform a maf file into vertexs[variant, allele]   """
import bmeg.enrichers.gene_enricher as gene_enricher
import logging

from bmeg.vertex import Callset, Gene
from bmeg.edge import AlleleCall

from bmeg.maf.maf_transform import main, get_value, MAFTransformer
from bmeg.maf.maf_transform import transform as parent_transform

CCLE_EXTENSION_CALLSET_KEYS = [
    'cDNA_Change', 'Codon_Change', 'Protein_Change',
    'isDeleterious', 'isTCGAhotspot', 'TCGAhsCnt',
    'isCOSMIChotspot', 'COSMIChsCnt', 'ExAC_AF', 'WES_AC',
    'WGS_AC', 'SangerWES_AC', 'SangerRecalibWES_AC', 'RNAseq_AC',
    'HC_AC', 'RD_AC'
]

CCLE_EXTENSION_MAF_KEYS = [
    'Genome_Change', 'Annotation_Transcript', 'cDNA_Change', 'Codon_Change', 'Protein_Change', 'isDeleterious', 'isTCGAhotspot', 'TCGAhsCnt',
    'isCOSMIChotspot', 'COSMIChsCnt', 'ExAC_AF', 'WES_AC', 'WGS_AC', 'SangerWES_AC', 'SangerRecalibWES_AC', 'RNAseq_AC', 'HC_AC', 'RD_AC',
]

TUMOR_SAMPLE_BARCODE = "Tumor_Sample_Barcode"  # 15
NORMAL_SAMPLE_BARCODE = "Matched_Norm_Sample_Barcode"  # 16


class CCLE_MAFTransformer(MAFTransformer):

    # callset source
    SOURCE = 'ccle'
    DEFAULT_PREFIX = SOURCE
    DEFAULT_MAF_FILE = 'source/ccle/CCLE_DepMap_18q3_maf_20180718.txt'

    def create_gene_gid(self, line):  # pragma nocover
        """ override, create gene_gid from line """
        symbol = line.get('Hugo_Symbol', None)
        try:
            genes = gene_enricher.get_gene(symbol)
            ensembl_id = genes[0]['ensembl_gene_id']
            return Gene.make_gid(gene_id=ensembl_id)
        except Exception as e:
            logging.warning(str(e))

    def callset_maker(self, allele, source, centerCol, method, line):
        """ create callset from line """
        aliquot_id = self.barcode_to_aliquot_id(line[TUMOR_SAMPLE_BARCODE])
        sample_callsets = []
        sample_calls = []
        callset = Callset(tumor_aliquot_id=aliquot_id,
                          normal_aliquot_id=None,
                          call_method=method,
                          source=source)
        sample_callsets.append(callset)
        sample_calls.append((self.allele_call_maker(allele, line),
                             callset.gid()))
        return sample_calls, sample_callsets

    def allele_call_maker(self, allele, line=None):
        """ create call from line """
        info = {}
        for k in CCLE_EXTENSION_CALLSET_KEYS:
            v = get_value(line, k, None)
            if v:
                info[k] = v
        return AlleleCall(info)

    def create_allele_dict(self, line, genome='GRCh37'):
        ''' return properly named allele dictionary, populated from line'''
        allele_dict = super(CCLE_MAFTransformer, self).create_allele_dict(line, genome)
        annotations = {}
        for key in CCLE_EXTENSION_MAF_KEYS:
            value = line.get(key, None)
            if value:
                annotations[key] = value

        allele_dict['annotations'].ccle = annotations
        return allele_dict


def transform(mafpath, prefix, emitter_name='json', skip=0):
    """ called from tests """
    return parent_transform(mafpath, prefix, CCLE_MAFTransformer.SOURCE, emitter_name, skip, transformer=CCLE_MAFTransformer())


if __name__ == '__main__':  # pragma: no cover
    main(transformer=CCLE_MAFTransformer())
