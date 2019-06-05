from glob import glob
import os
import pandas

import bmeg.ioutils
from bmeg.emitter import JSONEmitter
from bmeg import (Sample, Aliquot, Case, Project, Program,
                  Aliquot_Sample_Sample,
                  Sample_Case_Case,
                  Case_Projects_Project,
                  Sample_Projects_Project,
                  Aliquot_Projects_Project,
                  Case_Phenotypes_Phenotype,
                  Sample_Phenotypes_Phenotype,
                  Project_Programs_Program)
from bmeg.enrichers.phenotype_enricher import phenotype_factory


def transform(cellline_lookup_path="source/ccle/cellline_lookup.tsv",
              project_lookup_path="source/ccle/cellline_project_lookup.tsv",
              phenotype_lookup_path="source/ccle/cellline_phenotype_lookup.tsv",
              drug_response_path='source/ccle/CCLE_NP24.2009_Drug_data_2015.02.24.csv',
              expression_path="source/ccle/CCLE_depMap_19Q1_TPM.csv",
              maf_dir="source/ccle/mafs/*",
              emitter_prefix="ccle",
              emitter_directory="ccle"):

    celllines = bmeg.ioutils.read_lookup(cellline_lookup_path)
    projects = bmeg.ioutils.read_lookup(project_lookup_path)
    phenotypes = bmeg.ioutils.read_lookup(phenotype_lookup_path)

    emitter = JSONEmitter(directory=emitter_directory, prefix=emitter_prefix)
    prog = Program(submitter_id=Program.make_gid("CCLE"),
                   program_id="CCLE",
                   project_id="CCLE")
    emitter.emit_vertex(prog)

    raw_ids = {}
    drugs = {}

    input_stream = bmeg.ioutils.read_csv(drug_response_path)
    for line in input_stream:
        k = line['CCLE Cell Line Name']
        if k not in raw_ids:
            raw_ids[k] = None
        d = line['Compound']
        if k in drugs:
            drugs[k].append(d)
        else:
            drugs[k] = [d]

    input_stream = pandas.read_csv(expression_path, sep=",", index_col=0)
    for k, vals in input_stream.iterrows():
        if k not in raw_ids:
            raw_ids[k] = None

    for f in glob(maf_dir):
        k = os.path.basename(f).replace("_vs_NORMAL", "")
        if k not in raw_ids:
            raw_ids[k] = None

    emitted_celllines = {}
    emitted_projects = {}
    for i in raw_ids:
        emit_cellline = False
        if i in celllines:
            cellline_id = celllines[i]
        elif i.split("_")[0] in celllines:
            cellline_id = celllines[i.split("_")[0]]
        else:
            emit_cellline = True
            cellline_id = i

        if cellline_id in emitted_celllines:
            continue

        project_id = "CCLE_%s" % (projects.get(cellline_id, "Unknown"))
        proj = Project(submitter_id=Project.make_gid(project_id), project_id=project_id)
        if proj.gid() not in emitted_projects:
            emitter.emit_vertex(proj)
            emitter.emit_edge(
                Project_Programs_Program(
                    from_gid=proj.gid(),
                    to_gid=prog.gid(),
                ),
                emit_backref=True
            )
            emitted_projects[proj.gid()] = None

        c = Case(submitter_id=Case.make_gid(cellline_id),
                 case_id=cellline_id,
                 project_id=proj.gid())
        if emit_cellline:
            emitter.emit_vertex(c)
            # case <-> project edges
            emitter.emit_edge(
                Case_Projects_Project(
                    from_gid=c.gid(),
                    to_gid=proj.gid()
                ),
                emit_backref=True
            )

        sample_id = "CCLE:%s" % (cellline_id)
        s = Sample(submitter_id=Sample.make_gid(sample_id),
                   sample_id=sample_id,
                   project_id=proj.gid())
        emitter.emit_vertex(s)
        # sample <-> case edges
        emitter.emit_edge(
            Sample_Case_Case(
                from_gid=s.gid(),
                to_gid=c.gid()
            ),
            emit_backref=True
        )
        # sample <-> project edges
        emitter.emit_edge(
            Sample_Projects_Project(
                from_gid=s.gid(),
                to_gid=proj.gid()
            ),
            emit_backref=True
        )

        phenotype_name = phenotypes.get(cellline_id, None)
        if phenotype_name:
            pheno = phenotype_factory(phenotype_name)
            emitter.emit_vertex(pheno)
            # case <-> phenotype edges
            emitter.emit_edge(
                Case_Phenotypes_Phenotype(
                    from_gid=c.gid(),
                    to_gid=pheno.gid()
                ),
                emit_backref=True
            )
            # sample <-> phenotype edges
            emitter.emit_edge(
                Sample_Phenotypes_Phenotype(
                    from_gid=s.gid(),
                    to_gid=pheno.gid()
                ),
                emit_backref=True
            )

        def emit_aliquot(emitter, a, s, proj):
            # aliquot
            emitter.emit_vertex(a)
            # aliquot <-> sample edges
            emitter.emit_edge(
                Aliquot_Sample_Sample(
                    from_gid=a.gid(),
                    to_gid=s.gid()
                ),
                emit_backref=True
            )
            # aliquot <-> project edges
            emitter.emit_edge(
                Aliquot_Projects_Project(
                    from_gid=a.gid(),
                    to_gid=proj.gid()
                ),
                emit_backref=True
            )
            return

        for experiement_type in ["DrugResponse", "TranscriptExpression", "GeneExpression", "Callset"]:
            if experiement_type == "DrugResponse":
                for drug in drugs.get(i, []):
                    aliquot_id = "CCLE:%s:%s:%s" % (cellline_id, experiement_type, drug)
                    a = Aliquot(submitter_id=Aliquot.make_gid(aliquot_id),
                                aliquot_id=aliquot_id,
                                project_id=proj.gid())
                    emit_aliquot(emitter, a, s, proj)
            else:
                aliquot_id = "CCLE:%s:%s" % (cellline_id, experiement_type)
                a = Aliquot(submitter_id=Aliquot.make_gid(aliquot_id),
                            aliquot_id=aliquot_id,
                            project_id=proj.gid())
                emit_aliquot(emitter, a, s, proj)

        emitted_celllines[cellline_id] = None

    emitter.close()


if __name__ == '__main__':  # pragma: no cover
    transform()
