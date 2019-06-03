
import os
import pytest
import contextlib
from transform.ctrp.drug_response import transform
from bmeg.vertex import DrugResponse, Aliquot, Compound, Project


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
def cellline_lookup_path(request):
    return os.path.join(request.fspath.dirname, 'source/ccle/cellline_lookup.tsv')


@pytest.fixture
def project_lookup_path(request):
    return os.path.join(request.fspath.dirname, 'source/ccle/cellline_project_lookup.tsv')


def validate(helpers, emitter_directory, cellline_lookup_path, project_lookup_path, metadrugPath,
             metacelllinePath, responsePath, metaexperimentPath, curvePath):
    """ run xform and test results"""
    profile_file = os.path.join(emitter_directory, 'ctrp.DrugResponse.Vertex.json.gz')
    profile_in_file = os.path.join(emitter_directory, 'ctrp.ResponseIn.Edge.json.gz')
    response_to_file = os.path.join(emitter_directory, 'ctrp.ResponseTo.Edge.json.gz')
    compound_file = os.path.join(emitter_directory, 'ctrp.Compound.Vertex.json.gz')
    tested_in_file = os.path.join(emitter_directory, 'ctrp.TestedIn.Edge.json.gz')

    all_files = [profile_file, profile_in_file, response_to_file, compound_file, tested_in_file]

    # remove output
    with contextlib.suppress(FileNotFoundError):
        for f in all_files:
            os.remove(f)

    transform(
        cellline_lookup_path=cellline_lookup_path,
        project_lookup_path=project_lookup_path,
        metadrugPath=metadrugPath,
        metacelllinePath=metacelllinePath,
        responsePath=responsePath,
        metaexperimentPath=metaexperimentPath,
        curvePath=curvePath,
        emitter_directory=emitter_directory
    )

    # ratify
    helpers.assert_vertex_file_valid(DrugResponse, profile_file)
    helpers.assert_edge_file_valid(DrugResponse, Aliquot, profile_in_file)
    helpers.assert_edge_file_valid(DrugResponse, Compound, response_to_file)
    helpers.assert_vertex_file_valid(Compound, compound_file)
    helpers.assert_edge_file_valid(Compound, Project, tested_in_file)
    helpers.assert_edge_joins_valid(all_files, exclude_labels=['Aliquot', 'Project'])


def test_simple(helpers, emitter_directory, cellline_lookup_path, project_lookup_path, metadrugPath,
                metacelllinePath, responsePath, metaexperimentPath, curvePath):

    validate(helpers, emitter_directory, cellline_lookup_path, project_lookup_path, metadrugPath,
             metacelllinePath, responsePath, metaexperimentPath, curvePath)
