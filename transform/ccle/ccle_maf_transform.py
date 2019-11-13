import pandas

import bmeg.ioutils
from bmeg import Aliquot, Project, SomaticCallset, SomaticCallset_Alleles_Allele, SomaticCallset_Aliquots_Aliquot
from bmeg.emitter import new_emitter
from bmeg.maf.allele import make_minimal_allele, make_variant_call_data


def transform(mafpath="source/ccle/CCLE_DepMap_18q3_maf_20180718.txt",
              cellline_lookup_path="source/ccle/cellline_lookup.tsv",
              project_lookup_path="source/ccle/cellline_project_lookup.tsv",
              emitter_directory="ccle",
              emitter_prefix="maf"):

    id_lookup = bmeg.ioutils.read_lookup(cellline_lookup_path)
    project_lookup = bmeg.ioutils.read_lookup(project_lookup_path)

    emitter = new_emitter(directory=emitter_directory, prefix=emitter_prefix)

    emitted_alleles = {}
    emitted_callsets = {}
    for line in pandas.read_csv(mafpath, sep='\t', comment='#', dtype=str, chunksize=1):
        line = line.iloc[0, :].dropna().to_dict()

        allele = make_minimal_allele(line)
        if allele.gid() not in emitted_alleles:
            emitter.emit_vertex(allele)
            emitted_alleles[allele.gid()] = None

        tumor_aliquot_id = id_lookup.get(line.get('Broad_ID'), line.get('Broad_ID'))
        project_id = project_lookup.get(tumor_aliquot_id, None)

        callset = SomaticCallset(
            id=SomaticCallset.make_gid("CCLE", tumor_aliquot_id, None),
            tumor_aliquot_id=tumor_aliquot_id,
            normal_aliquot_id=None,
            project_id=Project.make_gid(project_id)
        )
        if callset.gid() not in emitted_callsets:
            emitter.emit_vertex(callset)
            emitted_callsets[callset.gid()] = None
            if callset.normal_aliquot_id:
                emitter.emit_edge(
                    SomaticCallset_Aliquots_Aliquot(
                        from_gid=callset.gid(),
                        to_gid=Aliquot.make_gid(callset.normal_aliquot_id),
                    ),
                    emit_backref=True
                )
            if callset.tumor_aliquot_id:
                emitter.emit_edge(
                    SomaticCallset_Aliquots_Aliquot(
                        from_gid=callset.gid(),
                        to_gid=Aliquot.make_gid(callset.tumor_aliquot_id),
                    ),
                    emit_backref=True
                )

        emitter.emit_edge(
            SomaticCallset_Alleles_Allele(
                from_gid=callset.gid(),
                to_gid=allele.gid(),
                data=make_variant_call_data(line, "Unknown")
            ),
            emit_backref=True
        )

    emitter.close()


if __name__ == '__main__':  # pragma: no cover
    transform()
