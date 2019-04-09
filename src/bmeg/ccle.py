import ujson
from types import SimpleNamespace as SN

import bmeg.ioutils
from bmeg.vertex import Aliquot, Sample, Case, Project
from bmeg.edge import HasAliquot, HasSample, HasCase
from bmeg.emitter import JSONEmitter, DebugEmitter


def build_project_lookup(ccle_sample_path="outputs/ccle/Sample.Vertex.json.gz"):
    projects = {}
    input_stream = bmeg.ioutils.reader(ccle_sample_path)
    for line in input_stream:
        sample = SN(**ujson.loads(line))
        ccle_attributes = sample.data['ccle_attributes']
        project_id = "_".join(ccle_attributes['Primary Disease'].split())
        aliquot_id = ccle_attributes['CCLE_Name']
        name_parts = aliquot_id.split('_')
        name_start = 1
        if len(name_parts) == 1:
            name_start = 0
        alias = '_'.join(name_parts[name_start:])
        projects[ccle_attributes['CCLE_Name']] = project_id
        projects[ccle_attributes['DepMap_ID']] = project_id
        projects[project_id] = project_id
        projects[alias] = project_id
    assert len(projects.keys()), 'No projects found, does {} exist?'.format(ccle_sample_path)
    return projects


def build_ccle2depmap_conversion_table(ccle_sample_path="outputs/ccle/Sample.Vertex.json.gz"):
    """create lookup table of CCLE_Name and Aliases to DepMap_ID"""
    samples = {}
    input_stream = bmeg.ioutils.reader(ccle_sample_path)
    for line in input_stream:
        sample = SN(**ujson.loads(line))
        ccle_attributes = SN(**sample.data['ccle_attributes'])

        samples[ccle_attributes.DepMap_ID] = ccle_attributes.DepMap_ID
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

    assert len(samples.keys()), 'Failed to create CCLE name to DepMap ID lookup, does {} exist?'.format(ccle_sample_path)
    return samples


def missing_ccle_cellline_factory(emitter, missing_ids,
                                  project_gids=[],
                                  project_prefix="", project_id="",
                                  project_lookup={}):
    """ create stub aliquot, sample, indiviudal, project and their edges """
    if not (isinstance(emitter, JSONEmitter) or isinstance(emitter, DebugEmitter)):
        raise TypeError("expected emitter to be an emitter")
    if isinstance(missing_ids, str):
        missing_ids = [missing_ids]
    elif isinstance(missing_ids, list):
        pass
    else:
        raise TypeError("expected missing_ids to be a list or string")
    if not isinstance(project_gids, list):
        raise TypeError("expected project_gids to be a list of strings")
    if not isinstance(project_prefix, str):
        raise TypeError("expected project_prefix to be a string")
    if not isinstance(project_id, str):
        raise TypeError("expected project_id to be a string")

    for aliquot_id in missing_ids:
        project = project_id
        if not project:
            project = aliquot_id
            alias = "_".join(project.split('_')[1:])
            if project in project_lookup:
                project = project_lookup[project]
            elif alias in project_lookup:
                project = project_lookup[alias]
            else:
                project = "Unknown"

        if project_prefix:
            project = "{}_{}".format(project_prefix, project)

        p = Project(project_id=project)
        if p.gid() not in project_gids:
            emitter.emit_vertex(p)
            project_gids.append(p.gid())

        c = Case(case_id=aliquot_id)
        emitter.emit_vertex(c)
        emitter.emit_edge(
            HasCase(),
            to_gid=c.gid(),
            from_gid=p.gid(),
        )

        s = Sample(sample_id=aliquot_id)
        emitter.emit_vertex(s)
        emitter.emit_edge(
            HasSample(),
            to_gid=s.gid(),
            from_gid=c.gid(),
        )

        a = Aliquot(aliquot_id=aliquot_id)
        emitter.emit_vertex(a)
        emitter.emit_edge(
            HasAliquot(),
            to_gid=a.gid(),
            from_gid=s.gid(),
        )

    return
