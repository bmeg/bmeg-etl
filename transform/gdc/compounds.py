import glob

from bmeg import (Case, Case_Compound_Compound)
from bmeg.emitter import JSONEmitter
from bmeg.enrichers.drug_enricher import compound_factory
from bmeg.ioutils import read_tsv


def transform(compounds="source/gdc/compounds/*.tsv",
              emitter_prefix=None,
              emitter_directory="gdc"):
    """ the only way to get drugs is to download files and parse them"""

    emitter = JSONEmitter(directory=emitter_directory, prefix=emitter_prefix)

    for path in glob.glob(compounds):
        tsv_in = read_tsv(path)
        c = 0
        dedupe = {}
        dedupe_compounds = {}
        for line in tsv_in:
            c += 1
            if c < 3:  # skip three header lines
                continue
            bcr_patient_uuid = line['bcr_patient_uuid']
            pharmaceutical_therapy_drug_name = line['pharmaceutical_therapy_drug_name']
            case_gid = Case.make_gid(bcr_patient_uuid)
            t = (bcr_patient_uuid, pharmaceutical_therapy_drug_name)
            # have we seen this combination before?
            if t not in dedupe:
                compound = compound_factory(name=pharmaceutical_therapy_drug_name)
                compound_gid = compound.gid()
                # have we seen this compound before?
                if compound_gid not in dedupe_compounds:
                    emitter.emit_vertex(compound)
                    dedupe_compounds[compound_gid] = True
                    emitter.emit_edge(
                        Case_Compound_Compound(
                            from_gid=case_gid,
                            to_gid=compound_gid
                        ),
                        emit_backref=True
                    )
                    # TODO emit edge from compound to project
                    dedupe[t] = True
    emitter.close()


if __name__ == '__main__':  # pragma: no cover
    transform()
