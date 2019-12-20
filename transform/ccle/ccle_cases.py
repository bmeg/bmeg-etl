import json
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


def transform(cellline_lookup_path="source/ccle/cellline_id_lookup.tsv",
              properties_lookup_path="source/ccle/cellline_properties_lookup.tsv",
              phenotype_lookup_path="source/ccle/cellline_phenotype_lookup.tsv",
              pharmacodb_cells_path="source/pharmacodb/cells.csv",
              pharmacodb_experiments_path="source/pharmacodb/experiments.csv",
              expression_path="source/ccle/CCLE_depMap_19Q1_TPM.csv",
              maf_path="source/ccle/CCLE_DepMap_18q3_maf_20180718.txt",
              emitter_prefix="ccle",
              emitter_directory="ccle"):

    emitter = JSONEmitter(directory=emitter_directory, prefix=emitter_prefix)

    celllines = bmeg.ioutils.read_lookup(cellline_lookup_path)
    properties = bmeg.ioutils.read_lookup(properties_lookup_path)
    phenotypes = bmeg.ioutils.read_lookup(phenotype_lookup_path)

    prog = Program(id=Program.make_gid("CCLE"),
                   program_id="CCLE")
    emitter.emit_vertex(prog)

    proj = Project(id=Project.make_gid("CCLE"),
                   project_id="CCLE")
    emitter.emit_vertex(proj)
    emitter.emit_edge(
        Project_Programs_Program(
            from_gid=proj.gid(),
            to_gid=prog.gid(),
        ),
        emit_backref=True
    )

    raw_ids = {}

    pharmacodb = pandas.merge(
        pandas.read_csv(pharmacodb_experiments_path)[["cell_id", "dataset_id"]],
        pandas.read_csv(pharmacodb_cells_path)[["cell_id", "cell_name"]],
        on="cell_id"
    ).drop_duplicates()
    pharmacodb_ccle = pharmacodb[pharmacodb.dataset_id == 1]

    for k in pharmacodb_ccle.cell_name.tolist():
        if k not in raw_ids:
            raw_ids[k] = None

    input_stream = pandas.read_csv(expression_path, sep=",", index_col=0)
    for k, vals in input_stream.iterrows():
        if k not in raw_ids:
            raw_ids[k] = None

    maf = pandas.read_csv(maf_path, sep="\t")[["Tumor_Sample_Barcode", "Broad_ID"]].drop_duplicates()
    for k in maf.Broad_ID.tolist():
        if k not in raw_ids:
            raw_ids[k] = None

    emitted_celllines = {}
    emitted_phenotypes = {}
    for raw_id in raw_ids:
        if raw_id in celllines:
            cellline_id = celllines[raw_id]
        elif raw_id.split("_")[0] in celllines:
            cellline_id = celllines[raw_id.split("_")[0]]
        else:
            cellline_id = raw_id

        if cellline_id in emitted_celllines:
            continue

        props = properties.get(cellline_id, None)
        if props:
            props = json.loads(props)

        case_id = "CCLE:%s" % (cellline_id)
        c = Case(id=Case.make_gid(case_id),
                 submitter_id=raw_id,
                 case_id=cellline_id,
                 cellline_attributes=props,
                 project_id=proj.gid())
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
        s = Sample(id=Sample.make_gid(sample_id),
                   submitter_id=raw_id,
                   sample_id=cellline_id,
                   cellline_attributes=props,
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

        aliquot_id = "CCLE:%s" % (cellline_id)
        a = Aliquot(id=Aliquot.make_gid(aliquot_id),
                    submitter_id=raw_id,
                    aliquot_id=cellline_id,
                    cellline_attributes=props,
                    project_id=proj.gid())
        emitter.emit_vertex(a)
        # aliquot <-> Sample edges
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

        phenotype_name = phenotypes.get(cellline_id, None)
        if phenotype_name:
            pheno = phenotype_factory(phenotype_name)
            if pheno.gid() not in emitted_phenotypes:
                emitter.emit_vertex(pheno)
                emitted_phenotypes[pheno.gid()] = None
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

        emitted_celllines[cellline_id] = None

    emitter.close()


if __name__ == '__main__':  # pragma: no cover
    transform()
