
""" transform a maf file into vertexs[variant, allele]   """

from glob import glob
from bmeg.vertex import File, Aliquot
from bmeg.edge import DerivedFrom
from bmeg.emitter import new_emitter
from bmeg.maf.maf_transform import MAFTransformer
from bmeg.ccle import build_ccle2depmap_conversion_table
import yaml
import os


TUMOR_SAMPLE_BARCODE = "Tumor_Sample_Barcode"  # 15
NORMAL_SAMPLE_BARCODE = "Matched_Norm_Sample_Barcode"  # 16

BIOSAMPLE_CONVERSION_TABLE = {}
MISSING_CELL_LINES = []


class CCLE_FileTransformer(MAFTransformer):

    # callset source
    SOURCE = 'ccle'
    DEFAULT_PREFIX = SOURCE

    def barcode_to_aliquot_id(self, barcode):
        """ create ccle sample barcode """
        global BIOSAMPLE_CONVERSION_TABLE
        global MISSING_CELL_LINES
        ccle_name_from_path = self.current_path.split('/')[-2]
        ccle_name_from_path = ccle_name_from_path.replace('_vs_NORMAL', '')
        if ccle_name_from_path in BIOSAMPLE_CONVERSION_TABLE:
            return BIOSAMPLE_CONVERSION_TABLE[ccle_name_from_path]
        elif ccle_name_from_path.split("_")[0] in BIOSAMPLE_CONVERSION_TABLE:
            return BIOSAMPLE_CONVERSION_TABLE[ccle_name_from_path.split("_")[0]]
        else:
            if ccle_name_from_path not in MISSING_CELL_LINES:
                MISSING_CELL_LINES.append(ccle_name_from_path)
            return ccle_name_from_path


def transform(mafpath, ccle_biosample_path, dvc_file, emitter_directory="ccle", emitter_prefix="maf"):
    """ entry point """

    # ensure that we have a lookup from CCLE native id to depmap id
    global BIOSAMPLE_CONVERSION_TABLE
    global MISSING_CELL_LINES
    BIOSAMPLE_CONVERSION_TABLE = build_ccle2depmap_conversion_table(ccle_biosample_path)
    emitter = new_emitter(name='json', directory=emitter_directory)
    transformer = CCLE_FileTransformer()
    md5s = {}
    with open(dvc_file, 'r') as stream:
        provenance = yaml.load(stream)
        for d in provenance['deps']:
            md5s[d['path']] = d['md5']

    for f in glob(mafpath):
        transformer.current_path = f
        aliquot_id = transformer.barcode_to_aliquot_id(None)
        md5sum = md5s.get(f, None)
        path = None
        if md5sum:
            path = 's3://bmeg/dvc/{}/{}'.format(md5sum[:2], md5sum[2:])
            file_size = os.path.getsize(f)
        file_dict = {'name': f, 'md5': md5sum, 'path': path, 'size': file_size}
        file = File(**file_dict)
        emitter.emit_vertex(file)
        emitter.emit_edge(DerivedFrom(), file.gid(), Aliquot.make_gid(aliquot_id))

    emitter.close()


if __name__ == '__main__':  # pragma: no cover
    transform(mafpath="source/ccle/mafs/*/vep.maf",
              ccle_biosample_path="outputs/ccle/Biosample.Vertex.json.gz",
              dvc_file="outputs.ccle.maf.dvc")
