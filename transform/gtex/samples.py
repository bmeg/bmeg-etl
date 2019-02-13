from bmeg.emitter import JSONEmitter
from bmeg.ioutils import read_tsv
from bmeg.vertex import Aliquot, Sample, Case, Project, Program
from bmeg.edge import AliquotFor, SampleFor, InProject, InProgram


def extract_case_id(sample_id):
    return "-".join(sample_id.split("-")[:2])


emitter = JSONEmitter("gtex", prefix="gtex")
samples = read_tsv("source/gtex/GTEx_v7_Annotations_SampleAttributesDS.txt")
cases = read_tsv("source/gtex/GTEx_v7_Annotations_SubjectPhenotypesDS.txt")

prog = Program("GTEx", gdc_attributes={})
emitter.emit_vertex(prog)

for row in cases:
    case_id = row["SUBJID"]
    i = Case(
        case_id,
        gdc_attributes={},
        gtex_attributes=row,
    )
    emitter.emit_vertex(i)

project_ids = []
for row in samples:
    sample_id = row["SAMPID"]
    case_id = extract_case_id(sample_id)

    p = Project(project_id="GTEx_{}".format(row["SMTS"]))
    if p.gid() not in project_ids:
        emitter.emit_vertex(p)
        emitter.emit_edge(
            InProgram(),
            p.gid(),
            prog.gid()
        )
    emitter.emit_edge(
        InProject(),
        Case.make_gid(case_id),
        p.gid(),
    )

    s = Sample(
        sample_id,
        gtex_attributes=row,
    )
    emitter.emit_vertex(s)
    emitter.emit_edge(
        SampleFor(),
        s.gid(),
        Case.make_gid(case_id),
    )

    a = Aliquot(sample_id)
    emitter.emit_vertex(a)
    emitter.emit_edge(
        AliquotFor(),
        a.gid(),
        s.gid(),
    )


emitter.close()
