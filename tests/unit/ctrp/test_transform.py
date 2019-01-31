
import os
import pytest
import contextlib
from transform.ctrp.transform import transform
from bmeg.vertex import DrugResponse, Aliquot, Compound, Biosample, Individual, Project


@pytest.fixture
def metadrugPath(request):
    return os.path.join(request.fspath.dirname, 'source/ctrp/v20.meta.per_compound.txt')


@pytest.fixture
def metacelllinePath(request):
    return os.path.join(request.fspath.dirname, 'source/ctrp/v20.meta.per_cell_line.txt')


@pytest.fixture
def responsePath(request):
    return os.path.join(request.fspath.dirname, 'source/ctrp/v20.data.curves_post_qc.txt')


@pytest.fixture
def metaexperimentPath(request):
    return os.path.join(request.fspath.dirname, 'source/ctrp/v20.meta.per_experiment.txt')


@pytest.fixture
def curvePath(request):
    return os.path.join(request.fspath.dirname, 'source/ctrp/v20.data.per_cpd_post_qc.txt')


@pytest.fixture
def biosample_path(request):
    return os.path.join(request.fspath.dirname, 'outputs/ccle/Biosample.Vertex.json.gz')


def validate(helpers, emitter_directory, biosample_path, metadrugPath,
             metacelllinePath, responsePath, metaexperimentPath, curvePath):
    """ run xform and test results"""
    profile_file = os.path.join(
        emitter_directory, 'ctrp.DrugResponse.Vertex.json.gz')
    profile_in_file = os.path.join(
        emitter_directory, 'ctrp.ResponseIn.Edge.json.gz')
    response_to_file = os.path.join(
        emitter_directory, 'ctrp.ResponseTo.Edge.json.gz')
    compound_file = os.path.join(
        emitter_directory, 'ctrp.Compound.Vertex.json.gz')

    all_files = [profile_file, profile_in_file,
                 response_to_file, compound_file]
    # remove output
    with contextlib.suppress(FileNotFoundError):
        for f in all_files:
            os.remove(f)

    transform(biosample_path=biosample_path,
              metadrugPath=metadrugPath,
              metacelllinePath=metacelllinePath,
              responsePath=responsePath,
              metaexperimentPath=metaexperimentPath,
              curvePath=curvePath,
              emitter_directory=emitter_directory)
    # ratify
    helpers.assert_vertex_file_valid(DrugResponse, profile_file)
    helpers.assert_edge_file_valid(DrugResponse, Aliquot, profile_in_file)
    helpers.assert_edge_file_valid(DrugResponse, Compound, response_to_file)
    helpers.assert_vertex_file_valid(Compound, compound_file)
    helpers.assert_edge_joins_valid(all_files, exclude_labels=['Aliquot'])
    # missing vertexes
    for f in ['ctrp.Aliquot.Vertex.json.gz', 'ctrp.Biosample.Vertex.json.gz',
              'ctrp.Individual.Vertex.json.gz', 'ctrp.Project.Vertex.json.gz']:
        v = eval(f.split('.')[1])
        f = os.path.join(emitter_directory, f)
        helpers.assert_vertex_file_valid(v, f)
    # missing edges
    for f, v1, v2 in [('ctrp.AliquotFor.Edge.json.gz', Aliquot, Biosample),
                      ('ctrp.BiosampleFor.Edge.json.gz',
                       Biosample, Individual),
                      ('ctrp.InProject.Edge.json.gz', Individual, Project)]:
        f = os.path.join(emitter_directory, f)
        helpers.assert_edge_file_valid(v1, v2, f)


def test_simple(helpers, emitter_directory, biosample_path, metadrugPath,
                metacelllinePath, responsePath, metaexperimentPath, curvePath):

    validate(helpers, emitter_directory, biosample_path, metadrugPath,
             metacelllinePath, responsePath, metaexperimentPath, curvePath)
