import os
import contextlib
import pytest
import shutil

from bmeg.emitter import JSONEmitter
from transform.pubmed.pubmed import parse_pubmed


@pytest.fixture
def pubmed_file(request):
    """ get the full path of the test fixture """
    return os.path.join(request.fspath.dirname, 'pubmed_test.xml')


def validate(helpers, emitter_directory, pubmed_file):
    publication_file = os.path.join(emitter_directory, 'Publication.Vertex.json.gz')

    # remove output
    with contextlib.suppress(FileNotFoundError):
        shutil.rmtree(emitter_directory)

    emitter = JSONEmitter(directory=emitter_directory)
    with open(pubmed_file) as handle:
        parse_pubmed(handle, emitter)

    helpers.assert_vertex_file_valid(publication_file)


def test_simple(helpers, emitter_directory, pubmed_file):
    validate(helpers, emitter_directory, pubmed_file)
