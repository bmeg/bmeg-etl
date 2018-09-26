from bmeg.emitter import JSONEmitter
from bmeg.ioutils import read_tsv
from bmeg.vertex import Aliquot, Biosample, Project, Individual
from bmeg.edge import AliquotFor, BiosampleFor, InProject


def extract_individual_id(sample_id):
    return "-".join(sample_id.split("-")[:2])


emitter = JSONEmitter("gtex", prefix="gtex")
samples = read_tsv("source/gtex/GTEx_v7_Annotations_SampleAttributesDS.txt")
individuals = read_tsv("source/gtex/GTEx_v7_Annotations_SubjectPhenotypesDS.txt")

project = Project("gtex", gdc_attributes={})
emitter.emit_vertex(project)

for row in individuals:
    individual_id = row["SUBJID"]
    i = Individual(
        individual_id,
        gdc_attributes={},
        gtex_attributes=row,
    )
    emitter.emit_vertex(i)
    emitter.emit_edge(
        InProject(),
        i.gid(),
        project.gid(),
    )

for row in samples:
    sample_id = row["SAMPID"]
    individual_id = extract_individual_id(sample_id)
    b = Biosample(
        sample_id,
        gtex_attributes=row,
    )
    emitter.emit_vertex(b)
    emitter.emit_edge(
        BiosampleFor(),
        b.gid(),
        Individual.make_gid(individual_id),
    )

    a = Aliquot(sample_id)
    emitter.emit_vertex(a)
    emitter.emit_edge(
        AliquotFor(),
        a.gid(),
        b.gid(),
    )

emitter.close()
