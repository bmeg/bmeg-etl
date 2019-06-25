
""" transform a maf file into vertexs[variant, allele]   """

from glob import glob
import bmeg.ioutils
from bmeg import Callset, Gene, Project
from bmeg.emitter import new_emitter
from bmeg.maf.maf_transform import get_value, MAFTransformer


CCLE_EXTENSION_CALLSET_INT_KEYS = {
    't_depth': 't_depth',
    't_ref_count': 't_ref_count',
    't_alt_count': 't_alt_count',
    'n_depth': 'n_depth',
    'n_ref_count': 'n_ref_count',
    'n_alt_count': 'n_alt_count'
}

CCLE_EXTENSION_CALLSET_KEYS = {
    'Reference_Allele': 'ref',
    'Tumor_Seq_Allele1': 'alt',
    'FILTER': 'filter',
    'Gene': 'ensembl_gene',
    'Feature': 'ensembl_transcript'
}

TUMOR_SAMPLE_BARCODE = "Tumor_Sample_Barcode"  # 15
NORMAL_SAMPLE_BARCODE = "Matched_Norm_Sample_Barcode"  # 16

SAMPLE_CONVERSION_TABLE = {}


class CCLE_MAFTransformer(MAFTransformer):

    # callset source
    TUMOR_ALLELE = 'Tumor_Seq_Allele2'

    def set_file(self, path):
        self.DEFAULT_MAF_FILE = path

    def create_gene_gid(self, line):  # pragma nocover
        ensembl_id = line.get('Gene', None)
        return Gene.make_gid(ensembl_id)

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
        return cellline_id

    def callset_maker(self, line, method):
        """ create callset from line """
        global PROJECT_CONVERSION_TABLE
        cellline_id = self.barcode_to_aliquot_id(line[TUMOR_SAMPLE_BARCODE])
        tumor_aliquot_id = "CCLE:%s:Callset" % (cellline_id)
        project_id = "CCLE_%s" % (PROJECT_CONVERSION_TABLE.get(cellline_id, "Unknown"))
        callset = Callset(
            submitter_id=Callset.make_gid("CCLE", tumor_aliquot_id, None),
            tumor_aliquot_id=tumor_aliquot_id,
            normal_aliquot_id=None,
            project_id=Project.make_gid(project_id)
        )
        allele_call = (self.allele_call_maker(line, method), callset.gid())
        return [allele_call], [callset]

    def allele_call_maker(self, line, method):
        """ create call from line """
        info = {
            "methods": [method],
        }
        for k, kn in CCLE_EXTENSION_CALLSET_KEYS.items():
            info[kn] = get_value(line, k, None)
        for k, kn in CCLE_EXTENSION_CALLSET_INT_KEYS.items():
            val = get_value(line, k, None)
            if val == "." or val == "" or val is None:
                info[kn] = None
            else:
                info[kn] = int(val)
        if info['filter'] is None:
            info['filter'] = 'PASS'
        return info


def transform(mafpath="source/ccle/mafs/*/vep.maf",
              cellline_lookup_path="source/ccle/cellline_lookup.tsv",
              project_lookup_path="source/ccle/cellline_project_lookup.tsv",
              emitter_directory="ccle",
              emitter_prefix="maf"):

    # ensure that we have a lookup from CCLE native id to depmap id
    global SAMPLE_CONVERSION_TABLE
    SAMPLE_CONVERSION_TABLE = bmeg.ioutils.read_lookup(cellline_lookup_path)
    # ensure that we have a lookup from CCLE native id to project id
    global PROJECT_CONVERSION_TABLE
    PROJECT_CONVERSION_TABLE = bmeg.ioutils.read_lookup(project_lookup_path)

    transformer = CCLE_MAFTransformer()
    emitter = new_emitter(name="json",
                          directory=emitter_directory,
                          prefix=emitter_prefix)

    for f in glob(mafpath):
        transformer.maf_convert(
            skip=1,
            emitter=emitter,
            mafpath=f
        )

    emitter.close()


if __name__ == '__main__':  # pragma: no cover
    transform()
