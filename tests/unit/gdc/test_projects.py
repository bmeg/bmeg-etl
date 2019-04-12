import os
import contextlib
import pytest

from transform.gdc.projects import transform
from bmeg.vertex import Project, Program
from bmeg.emitter import JSONEmitter


@pytest.fixture
def emitter_path_prefix(request):
    """ get the full path of the test output """
    return os.path.join(request.fspath.dirname, 'test/')


def validate(helpers, emitter_path_prefix):
    """ run xform and test results"""
    project_file = '{}Project.Vertex.json.gz'.format(emitter_path_prefix)
    program_file = '{}Program.Vertex.json.gz'.format(emitter_path_prefix)
    inprogram_file = '{}HasProject.Edge.json.gz'.format(emitter_path_prefix)
    all_files = [project_file, program_file, inprogram_file]

    # remove output
    with contextlib.suppress(FileNotFoundError):
        for f in all_files:
            os.remove(f)

    # create output
    emitter = JSONEmitter(emitter_path_prefix)
    transform(emitter)
    emitter.close()
    helpers.assert_vertex_file_valid(Project, project_file)
    helpers.assert_vertex_file_valid(Program, program_file)
    helpers.assert_edge_file_valid(Program, Project, inprogram_file)
    # validate vertex for all edges exist
    helpers.assert_edge_joins_valid(
        all_files,
        exclude_labels=[]
    )


def test_simple(helpers, emitter_path_prefix):
    validate(helpers, emitter_path_prefix)
