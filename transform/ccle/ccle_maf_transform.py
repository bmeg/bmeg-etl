
""" transform a maf file into vertexs[variant, allele]   """

from glob import glob
import bmeg.ioutils
from bmeg.vertex import Callset, Gene
from bmeg.edge import AlleleCall
from bmeg.emitter import new_emitter
from bmeg.maf.maf_transform import get_value, MAFTransformer


CCLE_EXTENSION_CALLSET_INT_KEYS = {
    # 't_depth' : 't_depth',
    # 't_ref_count' : 't_ref_count',
    # 't_alt_count' : 't_alt_count',
    # 'n_depth' : 'n_depth',
    # 'n_ref_count' : 'n_ref_count',
    # 'n_alt_count' : 'n_alt_count'
}

CCLE_EXTENSION_CALLSET_KEYS = {
    'Reference_Allele': 'ref',
    'Tumor_Seq_Allele1': 'alt',
    'FILTER': 'filter',
    'Gene': 'ensembl_gene',
    'Transcript_ID': 'ensembl_transcript'
}

TUMOR_SAMPLE_BARCODE = "Tumor_Sample_Barcode"  # 15
NORMAL_SAMPLE_BARCODE = "Matched_Norm_Sample_Barcode"  # 16

SAMPLE_CONVERSION_TABLE = {}


class CCLE_MAFTransformer(MAFTransformer):

    # callset source
    SOURCE = 'CCLE'
    DEFAULT_PREFIX = SOURCE
    TUMOR_ALLELE = 'Tumor_Seq_Allele2'

    def set_file(self, path):
        self.DEFAULT_MAF_FILE = path

    def create_gene_gid(self, line):  # pragma nocover
        ensembl_id = line.get('Gene', None)
        return Gene.make_gid(gene_id=ensembl_id)

    def barcode_to_aliquot_id(self, barcode):
        """ create ccle sample barcode """
        global SAMPLE_CONVERSION_TABLE
        ccle_name_from_path = self.current_path.split('/')[-2]
        ccle_name_from_path = ccle_name_from_path.replace('_vs_NORMAL', '')

        if ccle_name_from_path in SAMPLE_CONVERSION_TABLE:
            cellline_id = SAMPLE_CONVERSION_TABLE[ccle_name_from_path]
        elif ccle_name_from_path.split("_")[0] in SAMPLE_CONVERSION_TABLE:
            cellline_id = SAMPLE_CONVERSION_TABLE[ccle_name_from_path.split("_")[0]]
        else:
            cellline_id = ccle_name_from_path

        return "CCLE:%s:Callset" % (cellline_id)

    def callset_maker(self, allele, source, centerCol, method, line):
        """ create callset from line """
        aliquot_id = self.barcode_to_aliquot_id(line[TUMOR_SAMPLE_BARCODE])
        callset = Callset(tumor_aliquot_id=aliquot_id,
                          normal_aliquot_id=None,
                          source=source)
        sample_call = (self.allele_call_maker(allele, line, method), callset.gid())

        return [sample_call], [callset]

    def allele_call_maker(self, allele, line, method):
        """ create call from line """
        info = {
            "t_depth": 0,
            "n_depth": 0,
            "t_ref_count": 0,
            "t_alt_count": 0,
            "n_ref_count": 0,
            "n_alt_count": 0,
            "methods": ["ccle"],
        }
        for k, kn in CCLE_EXTENSION_CALLSET_KEYS.items():
            info[kn] = get_value(line, k, None)
        if info['filter'] is None:
            info['filter'] = 'PASS'
        # for k, kn in CCLE_EXTENSION_CALLSET_INT_KEYS.items():
        #     info[kn] = int(get_value(line, k, None))
        return AlleleCall(**info)


def transform(mafpath="source/ccle/mafs/*/vep.maf",
              cellline_lookup_path="source/ccle/cellline_lookup.tsv",
              emitter_directory="ccle",
              emitter_prefix="maf"):

    # ensure that we have a lookup from CCLE native id to depmap id
    global SAMPLE_CONVERSION_TABLE
    SAMPLE_CONVERSION_TABLE = bmeg.ioutils.read_lookup(cellline_lookup_path)

    transformer = CCLE_MAFTransformer()
    emitter = new_emitter(name="json",
                          directory=emitter_directory,
                          prefix=emitter_prefix)

    for f in glob(mafpath):
        transformer.maf_convert(
            skip=1,
            emitter=emitter,
            mafpath=f,
            source="CCLE")

    emitter.close()


if __name__ == '__main__':  # pragma: no cover
    transform()
