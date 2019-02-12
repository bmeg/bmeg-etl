
import os
import contextlib
import pytest
from transform.ccle.expression_tatlow import transform
from bmeg.vertex import TranscriptExpression, Aliquot, Biosample, Case, Project


@pytest.fixture
def ccle_tpm(request):
    """ get the full path of the test output """
    return os.path.join(request.fspath.dirname, 'source/ccle/expression/CCLE_tpm.tsv.gz')


@pytest.fixture
def biosample_path(request):
    """ get the full path of the test output """
    return os.path.join(request.fspath.dirname, 'outputs/ccle/Biosample.Vertex.json.gz')


def validate(helpers, ccle_tpm, biosample_path, emitter_directory):
    """ run xform and test results"""
    expression_file = os.path.join(emitter_directory, 'tatlow.TranscriptExpression.Vertex.json.gz')
    expression_of_file = os.path.join(emitter_directory, 'tatlow.TranscriptExpressionOf.Edge.json.gz')
    all_files = [expression_file, expression_of_file]
    # remove output
    with contextlib.suppress(FileNotFoundError):
        for f in all_files:
            os.remove(f)

    # create output
    transform(source_path=ccle_tpm, biosample_path=biosample_path, emitter_directory=emitter_directory)
    # ratify
    helpers.assert_vertex_file_valid(TranscriptExpression, expression_file)
    helpers.assert_edge_file_valid(TranscriptExpression, Aliquot, expression_of_file)
    helpers.assert_edge_joins_valid(all_files, exclude_labels=['Aliquot'])

    # test extra cell lines
    helpers.assert_vertex_file_valid(Aliquot, os.path.join(emitter_directory, 'tatlow.Aliquot.Vertex.json.gz'))
    helpers.assert_vertex_file_valid(Case, os.path.join(emitter_directory, 'tatlow.Case.Vertex.json.gz'))
    helpers.assert_vertex_file_valid(Project, os.path.join(emitter_directory, 'tatlow.Project.Vertex.json.gz'))

    helpers.assert_edge_file_valid(Aliquot, Biosample, os.path.join(emitter_directory, 'tatlow.AliquotFor.Edge.json.gz'))
    helpers.assert_edge_file_valid(Biosample, Case, os.path.join(emitter_directory, 'tatlow.BiosampleFor.Edge.json.gz'))
    helpers.assert_edge_file_valid(Case, Project, os.path.join(emitter_directory, 'tatlow.InProject.Edge.json.gz'))


def test_simple(helpers, ccle_tpm, biosample_path, emitter_directory):
    """ just run validate"""
    validate(helpers, ccle_tpm, biosample_path, emitter_directory)
