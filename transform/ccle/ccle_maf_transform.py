"""
Transform a maf file into AlleleCall and Callset vertices.
"""

import bmeg.maf.gene_enricher as gene_enricher
import logging

from bmeg.vertex import Biosample, Callset, Gene
from bmeg.edge import AlleleCall
from bmeg.util.cli import default_argument_parser
from bmeg.util.logging import default_logging
from bmeg.emitter import JSONEmitter
from bmeg.maf.maf_transform import get_value, MAFTransformer


TUMOR_SAMPLE_BARCODE = "Tumor_Sample_Barcode"  # 15
NORMAL_SAMPLE_BARCODE = "Matched_Norm_Sample_Barcode"  # 16


class CCLE_MAFTransformer(MAFTransformer):

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
        sample = self.barcode_to_sampleid(line[TUMOR_SAMPLE_BARCODE])
        sample = Biosample.make_gid(sample)
        sample_callsets = []
        sample_calls = []
        callset = Callset(sample,
                          None,
                          method, source)
        sample_callsets.append(callset)
        sample_calls.append((self.allele_call_maker(allele, line),
                             callset.gid()))
        return sample_calls, sample_callsets

    def allele_call_maker(self, allele, line=None):
        """ create call from line """
        keys = ['cDNA_Change', 'Codon_Change', 'Protein_Change',
                'isDeleterious', 'isTCGAhotspot', 'TCGAhsCnt',
                'isCOSMIChotspot', 'COSMIChsCnt', 'ExAC_AF', 'WES_AC',
                'WGS_AC', 'SangerWES_AC', 'SangerRecalibWES_AC', 'RNAseq_AC',
                'HC_AC', 'RD_AC']
        info = {}
        for k in keys:
            v = get_value(line, k, None)
            if v:
                info[k] = v
        return AlleleCall(info)


parser = default_argument_parser()
options = parser.parse_args()
default_logging(options.loglevel)

emitter = JSONEmitter("ccle")

maf_path = "source/ccle/ccle2maf_081117.txt"
transformer = CCLE_MAFTransformer()
transformer.maf_convert(
    emitter=emitter,
    mafpath=maf_path,
    workers=10,
    source="ccle",
    gz=False,
)

emitter.close()
