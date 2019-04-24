from glob import glob
import os
import pandas

import bmeg.ioutils
from bmeg.emitter import JSONEmitter
from bmeg.vertex import Sample, Aliquot, Case, Project, Program
from bmeg.edge import AliquotFor, SampleFor, InProject, InProgram, PhenotypeOf
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
    prog = Program(program_id="CCLE")
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
        project_id = "CCLE_%s" % (projects.get(i, "Unknown"))
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

        p = Project(project_id=project_id)
        if p.gid() not in emitted_projects:
            emitter.emit_vertex(p)
            emitter.emit_edge(
                InProgram(),
                p.gid(),
                prog.gid(),
            )
            emitted_projects[p.gid()] = None

        c = Case(case_id=cellline_id)
        if emit_cellline:
            emitter.emit_vertex(c)
        emitter.emit_edge(
            InProject(),
            c.gid(),
            p.gid(),
        )

        s = Sample(sample_id=cellline_id)
        emitter.emit_vertex(s)
        emitter.emit_edge(
            SampleFor(),
            s.gid(),
            c.gid(),
        )
        emitter.emit_edge(
            InProject(),
            s.gid(),
            p.gid(),
        )

        phenotype_name = phenotypes.get(cellline_id, None)
        if phenotype_name:
            pheno = phenotype_factory(phenotype_name)
            emitter.emit_vertex(pheno)
            emitter.emit_edge(
                PhenotypeOf(),
                c.gid(),
                pheno.gid()
            )
            emitter.emit_edge(
                PhenotypeOf(),
                s.gid(),
                pheno.gid()
            )

        for experiement_type in ["DrugResponse", "TranscriptExpression", "GeneExpression", "Callset"]:            
            if experiement_type == "DrugResponse":
                for drug in drugs.get(i, []):
                    aliquot_id = "%s:%s:%s" % (cellline_id, experiement_type, drug)
                    a = Aliquot(aliquot_id=aliquot_id)
                    emitter.emit_vertex(a)
                    emitter.emit_edge(
                        AliquotFor(),
                        a.gid(),
                        s.gid(),
                    )
                    emitter.emit_edge(
                        InProject(),
                        a.gid(),
                        p.gid(),
                    )
            else:
                aliquot_id = "%s:%s" % (cellline_id, experiement_type)
                a = Aliquot(aliquot_id=aliquot_id)
                emitter.emit_vertex(a)
                emitter.emit_edge(
                    AliquotFor(),
                    a.gid(),
                    s.gid(),
                )
                emitter.emit_edge(
                    InProject(),
                    a.gid(),
                    p.gid(),
                )

        emitted_celllines[cellline_id] = None


if __name__ == '__main__':  # pragma: no cover
    transform()
