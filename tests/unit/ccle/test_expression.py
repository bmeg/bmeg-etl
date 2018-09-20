
import os
import contextlib
import pytest
from transform.ccle.expression import transform
from bmeg.vertex import Expression, Aliquot
import json


@pytest.fixture
def gct_file(request):
    """ get the full path of the test output """
    return os.path.join(request.fspath.dirname, 'source/ccle/CCLE_DepMap_18q3_RNAseq_RPKM_20180718.gct')


def validate(helpers, gct_file, emitter_directory):
    """ run xform and test results"""
    expression_file = os.path.join(emitter_directory, 'ccle.Expression.Vertex.json')
    expression_of_file = os.path.join(emitter_directory, 'ccle.ExpressionOf.Edge.json')

    all_files = [expression_file, expression_of_file]
    # remove output
    with contextlib.suppress(FileNotFoundError):
        for f in all_files:
            os.remove(f)

    # create output
    transform(path=gct_file, emitter_directory=emitter_directory)
    # ratify
    helpers.assert_vertex_file_valid(Expression, expression_file)
    helpers.assert_edge_file_valid(Expression, Aliquot, expression_of_file)
    helpers.assert_edge_joins_valid(all_files, exclude_labels=['Aliquot'])
    # ensure broad ids not used for aliquot
    with open(expression_of_file, 'r', encoding='utf-8') as f:
        for line in f:
            obj = json.loads(line)
            assert len(obj['to'].split()) == 1, 'Broad ids should not be in Aliquot id'


def test_simple(helpers, gct_file, emitter_directory):
    """ just run validate"""
    validate(helpers, gct_file, emitter_directory)
