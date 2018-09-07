
import os
import contextlib
import pytest
import json
from transform.gdsc.response import transform
from bmeg.vertex import Compound, DrugResponse, Biosample


@pytest.fixture
def emitter_path_prefix(request):
    """ get the full path of the test output """
    return os.path.join(request.fspath.dirname, 'test/test')


@pytest.fixture
def GDSC_AUC_file(request):
    """ get the full path of the test fixture """
    return os.path.join(request.fspath.dirname, 'source/gdsc/GDSC_AUC.csv')


ALL_FILES = """
Compound.Vertex.json
DrugResponse.Vertex.json
DrugResponseIn.Edge.json
ResponseTo.Edge.json
""".strip().split()


def validate(helpers, emitter_path_prefix, GDSC_AUC_file):

    all_files = ['{}.{}'.format(emitter_path_prefix, f) for f in ALL_FILES]
    # remove output
    with contextlib.suppress(FileNotFoundError):
        for f in all_files:
            os.remove(f)

    # create output
    transform(path=GDSC_AUC_file, prefix=emitter_path_prefix)

    compounds = all_files[0]
    drug_responses = all_files[1]
    drug_response_ins = all_files[2]
    response_tos = all_files[3]

    helpers.assert_vertex_file_valid(Compound, compounds)
    helpers.assert_vertex_file_valid(DrugResponse, drug_responses)
    helpers.assert_edge_file_valid(DrugResponse, Biosample, drug_response_ins)
    helpers.assert_edge_file_valid(DrugResponse, Compound, response_tos)

    # validate vertex for all edges exist
    helpers.assert_edge_joins_valid(all_files, exclude_labels=['Biosample'])

    # test Compound contents
    c = 0
    with open(compounds, 'r', encoding='utf-8') as f:
        for line in f:
            json.loads(line)
            c += 1
    assert c == 9, 'Should have 9 compounds'


def test_simple(helpers, emitter_path_prefix, GDSC_AUC_file):
    validate(helpers, emitter_path_prefix, GDSC_AUC_file)
