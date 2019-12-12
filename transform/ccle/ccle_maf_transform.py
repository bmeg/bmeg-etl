import pandas

import bmeg.ioutils
from bmeg import Aliquot, Project, SomaticCallset, SomaticCallset_Alleles_Allele, SomaticCallset_Aliquots_Aliquot
from bmeg.emitter import new_emitter
from bmeg.maf import make_minimal_allele, make_variant_call_data


def transform(mafpath="source/ccle/CCLE_DepMap_18q3_maf_20180718.txt",
              cellline_lookup_path="source/ccle/cellline_id_lookup.tsv",
              emitter_directory="ccle",
              emitter_prefix="maf"):

    id_lookup = bmeg.ioutils.read_lookup(cellline_lookup_path)

    emitter = new_emitter(directory=emitter_directory, prefix=emitter_prefix)

    emitted_alleles = {}
    emitted_callsets = {}
    maf = pandas.read_csv(mafpath, sep='\t', comment='#', dtype=str)
    for index, line in maf.iterrows():
        line = line.dropna().to_dict()

        allele = make_minimal_allele(line, alternate_bases='Tumor_Seq_Allele1')
        if allele.gid() not in emitted_alleles:
            emitter.emit_vertex(allele)
            emitted_alleles[allele.gid()] = None

        tumor_aliquot_id = id_lookup.get(line.get('Broad_ID'), line.get('Broad_ID'))

        callset = SomaticCallset(
            id=SomaticCallset.make_gid("CCLE", tumor_aliquot_id, None),
            tumor_aliquot_id=tumor_aliquot_id,
            normal_aliquot_id=None,
            project_id=Project.make_gid("CCLE")
        )
        if callset.gid() not in emitted_callsets:
            emitter.emit_vertex(callset)
            emitted_callsets[callset.gid()] = None
            emitter.emit_edge(
                SomaticCallset_Aliquots_Aliquot(
                    from_gid=callset.gid(),
                    to_gid=Aliquot.make_gid("CCLE:%s" % (callset.tumor_aliquot_id)),
                ),
                emit_backref=True
            )

        emitter.emit_edge(
            SomaticCallset_Alleles_Allele(
                from_gid=callset.gid(),
                to_gid=allele.gid(),
                data=make_variant_call_data(line, alternate_bases='Tumor_Seq_Allele1')
            ),
            emit_backref=True
        )

    emitter.close()


if __name__ == '__main__':  # pragma: no cover
    transform()
