import ujson
from types import SimpleNamespace as SN

import bmeg.ioutils
from bmeg.vertex import Aliquot, Biosample, Individual, Project
from bmeg.edge import AliquotFor, BiosampleFor, InProject
from bmeg.emitter import JSONEmitter, DebugEmitter


def build_ccle2depmap_conversion_table(ccle_biosample_path="outputs/ccle/Biosample.Vertex.json.gz"):
    """create lookup table of CCLE_Name and Aliases to DepMap_ID"""
    samples = {}
    input_stream = bmeg.ioutils.reader(ccle_biosample_path)
    for line in input_stream:
        biosample = SN(**ujson.loads(line))
        ccle_attributes = SN(**biosample.data['ccle_attributes'])
        samples[ccle_attributes.CCLE_Name] = ccle_attributes.DepMap_ID
        if "NORMAL" in ccle_attributes.CCLE_Name or "MERGED" in ccle_attributes.CCLE_Name:
            continue
        # add prefix (e.g. CAL62 instead of CAL62_THYROID)
        prefix = ccle_attributes.CCLE_Name.split('_')[0]
        aliases = [prefix]
        # aliases = [prefix, ccle_attributes.Aliases]
        # split_aliases = ccle_attributes.Aliases.split(";")
        # if len(split_aliases) > 1:
        #     for x in split_aliases:
        #         aliases.append(x)
        # exclude non uniq aliases
        blacklisted_aliases = ["KMH2", "MS1", "TT"]
        for a in aliases:
            if a not in blacklisted_aliases:
                samples[a] = ccle_attributes.DepMap_ID

    assert len(samples.keys()), 'No biosamples found, does {} exist?'.format(ccle_biosample_path)
    return samples


def missing_ccle_cellline_factory(emitter, source, missing_ids, project_id=""):
    """ create stub aliquot, biosample, indiviudal, project and their edges """
    if not (isinstance(emitter, JSONEmitter) or isinstance(emitter, DebugEmitter)):
        raise TypeError("expected emitter to be an emitter")
    if not isinstance(source, str):
        raise TypeError("expected source to be a string")
    if not isinstance(project_id, str):
        raise TypeError("expected project_id to be a string")
    if isinstance(missing_ids, str):
        missing_ids = [missing_ids]
    elif isinstance(missing_ids, list):
        pass
    else:
        raise TypeError("expected missing_ids to be a list or string")

    individual_gids = project_gids = []
    for aliquot_id in missing_ids:
        project = project_id
        if not project:
            project = aliquot_id
            # strip off prefix
            name_parts = project.split('_')
            name_start = 1
            if len(name_parts) == 1:
                name_start = 0
            project = '{}:{}'.format(source, '_'.join(name_parts[name_start:]))

        aliquot_id = '{}:{}'.format(source, aliquot_id)

        p = Project(project_id=project)
        if p.gid() not in project_gids:
            emitter.emit_vertex(p)
            project_gids.append(p.gid())

        i = Individual(individual_id=aliquot_id)
        if i.gid() not in individual_gids:
            emitter.emit_vertex(i)
            emitter.emit_edge(
                InProject(),
                i.gid(),
                p.gid(),
            )
            individual_gids.append(i.gid())

        b = Biosample(biosample_id=aliquot_id)
        emitter.emit_vertex(b)
        emitter.emit_edge(
            BiosampleFor(),
            b.gid(),
            i.gid(),
        )

        a = Aliquot(aliquot_id=aliquot_id)
        emitter.emit_vertex(a)
        emitter.emit_edge(
            AliquotFor(),
            a.gid(),
            b.gid(),
        )

    return
