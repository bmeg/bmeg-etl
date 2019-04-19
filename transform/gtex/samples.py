from bmeg.emitter import JSONEmitter
from bmeg.ioutils import read_tsv
from bmeg import Aliquot, Sample, Case, Project, Program
#from bmeg.edge import HasAliquot, HasSample, HasCase, HasProject


def extract_case_id(sample_id):
    return "-".join(sample_id.split("-")[:2])


# https://github.com/guigolab/DeepFryer/blob/master/Description/sample.tsv
def sample_meta_map(meta):
    out = {
        "biospecimen_anatomic_site" : meta["SMTS"]
    }
    return out


emitter = JSONEmitter("gtex", prefix="gtex")
samples = read_tsv("source/gtex/GTEx_v7_Annotations_SampleAttributesDS.txt")
cases = read_tsv("source/gtex/GTEx_v7_Annotations_SubjectPhenotypesDS.txt")

prog = Program(id="GTEx")
emitter.emit_vertex(prog)

case_data = {}
for row in cases:
    case_id = row["SUBJID"]
    case_data[case_id] = row


case_ids = []
project_ids = []
for row in samples:
    sample_id = row["SAMPID"]
    case_id = extract_case_id(sample_id)
    project_id = "GTEx_{}".format(row["SMTS"])
    p = Project(id=project_id, code="", dbgap_accession_number="", name=row["SMTS"], program_id="GTEx")
    if p.id not in project_ids:
        emitter.emit_vertex(p)
        #emitter.emit_edge(
        #    HasProject(),
        #    to_gid=p.gid(),
        #    from_gid=prog.gid()
        #)
    #emitter.emit_edge(
    #    HasCase(),
    #    to_gid=Case.make_gid(case_id),
    #    from_gid=p.gid(),
    #)

    if case_id not in case_ids:
        i = Case(
            id=case_id,
            project_id=p.id,
            **case_data[case_id],
        )
        emitter.emit_vertex(i)
        case_ids.append(case_id)
    
    s = Sample(
        id=sample_id,
        cases=[case_id],
        submitter_id=sample_id,
        type="",
        **sample_meta_map(row)
    )
    emitter.emit_vertex(s)
    #emitter.emit_edge(
    #    HasSample(),
    #    to_gid=s.gid(),
    #    from_gid=Case.make_gid(case_id),
    #)

    a = Aliquot(id=sample_id, type="", submitter_id=sample_id)
    emitter.emit_vertex(a)
    #emitter.emit_edge(
    #    HasAliquot(),
    #    to_gid=a.gid(),
    #    from_gid=s.gid(),
    #)


emitter.close()
