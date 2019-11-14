from glob import glob

import bmeg.ioutils
from bmeg import Aliquot, Project, SomaticCallset, SomaticCallset_Alleles_Allele, SomaticCallset_Aliquots_Aliquot
from bmeg.emitter import new_emitter
from bmeg.vcf import make_minimal_allele, make_variant_call_data


def transform(vcf_dir="source/cellmodelpassports/vcfs/*",
              cellline_lookup_path="source/ccle/cellline_id_lookup.tsv",
              emitter_prefix=None,
              emitter_directory="gdsc"):

    id_lookup = bmeg.ioutils.read_lookup(cellline_lookup_path)

    emitter = new_emitter(directory=emitter_directory, prefix=emitter_prefix)

    vcf_files = glob(vcf_dir)
    vcf_files.sort()

    emitted_alleles = {}
    emitted_callsets = {}
    for f in vcf_files:
        pass

    emitter.close()


if __name__ == '__main__':  # pragma: no cover
    transform()
