import logging
import pandas

import bmeg.ioutils
from bmeg import Aliquot, Project, SomaticCallset, SomaticCallset_Alleles_Allele, SomaticCallset_Aliquots_Aliquot
from bmeg.emitter import new_emitter
from bmeg.maf import make_minimal_allele, make_variant_call_data


def transform(mafpath='source/mc3/mc3.v0.2.8.PUBLIC.maf.gz',
              id_lookup_path='source/gdc/id_lookup.tsv',
              project_lookup_path='source/gdc/project_lookup.tsv',
              emitter_directory='mc3',
              emitter_prefix=None):

    id_lookup = bmeg.ioutils.read_lookup(id_lookup_path)
    project_lookup = bmeg.ioutils.read_lookup(project_lookup_path)

    emitter = new_emitter(directory=emitter_directory, prefix=emitter_prefix)

    emitted_alleles = {}
    emitted_callsets = {}
    maf = pandas.read_csv(mafpath, sep='\t', comment='#', dtype=str, low_memory=False,
                          na_values=['.', 'null', 'NA', 'N/A'],
                          keep_default_na=False)
    for index, row in maf.iterrows():
        line = row.dropna().to_dict()

        allele = make_minimal_allele(line)
        if allele.gid() not in emitted_alleles:
            emitter.emit_vertex(allele)
            emitted_alleles[allele.gid()] = None

        tumor_aliquot_gid = id_lookup.get(line.get('Tumor_Sample_Barcode'), line.get('Tumor_Sample_Barcode'))
        normal_aliquot_gid = id_lookup.get(line.get('Matched_Norm_Sample_Barcode'), line.get('Matched_Norm_Sample_Barcode'))
        project_id = project_lookup.get(tumor_aliquot_gid, None)
        callset = SomaticCallset(
            id=SomaticCallset.make_gid("MC3", tumor_aliquot_gid, normal_aliquot_gid),
            tumor_aliquot_id=tumor_aliquot_gid,
            normal_aliquot_id=normal_aliquot_gid,
            project_id=Project.make_gid(project_id)
        )
        if callset.gid() not in emitted_callsets:
            emitter.emit_vertex(callset)
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
            emitted_callsets[callset.gid()] = None

        call_methods = line.get('CENTERS', '')
        call_methods = [call_method.replace('*', '') for call_method in call_methods.split("|")]
        emitter.emit_edge(
            SomaticCallset_Alleles_Allele(
                from_gid=callset.gid(),
                to_gid=allele.gid(),
                data=make_variant_call_data(line, call_methods)
            ),
            emit_backref=True
        )

    emitter.close()


if __name__ == '__main__':  # pragma: no cover
    transform()
