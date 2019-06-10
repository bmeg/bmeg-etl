import bmeg.ioutils
from bmeg.emitter import JSONEmitter
from bmeg import Case


def transform(path="source/ccle/DepMap-2019q1-celllines.csv_v2.csv",
              emitter_prefix="depmap",
              emitter_directory="ccle"):

    emitter = JSONEmitter(directory=emitter_directory, prefix=emitter_prefix)
    reader = bmeg.ioutils.read_csv(path)

    # Load sample metadata.
    case_gids = {}
    # DepMap_ID,CCLE_Name,Aliases,COSMIC_ID,Sanger ID,Primary Disease,Subtype Disease,Gender,Source
    # ACH-000001,NIHOVCAR3_OVARY,NIH:OVCAR-3;OVCAR3,905933,2201,Ovarian Cancer,"Adenocarcinoma, high grade serous",Female,ATCC
    for row in reader:
        if "MERGED" in row["CCLE_Name"]:
            continue

        props = {}
        for k, v in row.items():
            props["_".join(k.split())] = v

        c = Case(
            submitter_id=Case.make_gid(row["DepMap_ID"]),
            case_id=row["DepMap_ID"],
            cellline_attributes=props,
            project_id=''
        )
        if c.gid() not in case_gids:
            emitter.emit_vertex(c)
            case_gids[c.gid()] = None

    emitter.close()


if __name__ == '__main__':  # pragma: no cover
    transform()
