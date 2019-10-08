import glob
import re

from bmeg import (Case, Project, Case_Compounds_Compound, Compound_Projects_Project)
from bmeg.emitter import JSONEmitter
from bmeg.enrichers.drug_enricher import compound_factory
from bmeg.ioutils import read_tsv, read_lookup


def transform(compounds="source/gdc/compounds/*.tsv",
              id_lookup_path="source/gdc/id_lookup.tsv",
              project_lookup_path="source/gdc/project_lookup.tsv",
              emitter_prefix="gdc",
              emitter_directory="gdc"):

    ids = read_lookup(id_lookup_path)
    projects = read_lookup(project_lookup_path)

    emitter = JSONEmitter(directory=emitter_directory, prefix=emitter_prefix)

    dedup_proj = {}
    dedup_case = {}
    dedup_compounds = {}
    for path in glob.glob(compounds):
        tsv_in = read_tsv(path)
        c = 0
        for line in tsv_in:
            c += 1
            if c < 3:  # skip three header lines
                continue
            bcr_patient_barcode = line["bcr_patient_barcode"]
            drug_name = line["pharmaceutical_therapy_drug_name"].lower()
            case_gid = Case.make_gid(ids.get(bcr_patient_barcode, bcr_patient_barcode))
            project_gid = Project.make_gid(projects.get(bcr_patient_barcode, None))

            if any(["unknown" in drug_name,
                    drug_name == "yes",
                    drug_name == "no",
                    drug_name == "[not available]",
                    drug_name == "not specified",
                    drug_name == "not otherwise specified"]):
                continue

            drug_name = re.search("([A-Za-z0-9-_ ]+)(\(.*\))?", drug_name).group(1).strip()
            cpd_names = [x.strip() for x in re.split(",|\+", drug_name)]
            for cpd_name in cpd_names:
                if "placebo" in cpd_name:
                    continue
                compound = compound_factory(name=cpd_name)
                compound_gid = compound.gid()

                # have we seen this compound before?
                if compound_gid not in dedup_compounds:
                    emitter.emit_vertex(compound)
                    dedup_compounds[compound_gid] = True

                # have we seen this combination before?
                t = (case_gid, compound_gid)
                if t not in dedup_case:
                    emitter.emit_edge(
                        Case_Compounds_Compound(
                            from_gid=case_gid,
                            to_gid=compound_gid
                        ),
                        emit_backref=True
                    )
                    dedup_case[t] = True

                # have we seen this combination before?
                t = (project_gid, compound_gid)
                if t not in dedup_proj:
                    emitter.emit_edge(
                        Compound_Projects_Project(
                            from_gid=compound_gid,
                            to_gid=project_gid
                        ),
                        emit_backref=True
                    )
                    dedup_proj[t] = True

    emitter.close()


if __name__ == "__main__":  # pragma: no cover
    transform()
