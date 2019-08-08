import contextlib
import os
import shutil
import pytest
from transform.celllines.link_cases import transform


@pytest.fixture
def case_files(request):
    return [os.path.join(request.fspath.dirname, 'outputs/ccle/ccle.Case.Vertex.json'),
            os.path.join(request.fspath.dirname, 'outputs/ctrp/ctrp.Case.Vertex.json'),
            os.path.join(request.fspath.dirname, 'outputs/gdsc/gdsc.Case.Vertex.json')]


def validate(helpers, emitter_directory, case_files):
    """ run xform and test results"""

    edge_file = os.path.join(emitter_directory, 'Case_SameAs_Case.Edge.json.gz')

    # remove output
    with contextlib.suppress(FileNotFoundError):
        shutil.rmtree(emitter_directory)

    # create output
    transform(
        case_files=case_files,
        emitter_prefix=None,
        emitter_directory=emitter_directory
    )

    n = helpers.assert_edge_file_valid(edge_file)
    assert n == 34, "expected 34 edges"

    # all_files = case_files + [edge_file]
    # helpers.assert_edge_joins_valid(
    #     all_files,
    #     exclude_labels=[]
    # )


def test_simple(helpers, emitter_directory, case_files):
    validate(helpers, emitter_directory, case_files)
