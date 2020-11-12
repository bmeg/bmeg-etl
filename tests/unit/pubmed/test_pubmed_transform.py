import os
import contextlib
import pytest
import shutil

from bmeg import Publication, Project

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

    outputs = []
    with open(pubmed_file) as handle:
        parse_pubmed(handle, outputs)

    emitter = JSONEmitter(directory=emitter_directory)
    for o in outputs:
        p = Publication(id=Publication.make_gid(o["url"]),
                        url=o["url"], title=o["title"], abstract=o["abstract"],
                        text="", date=o["date"], author=o["author"], citation=[],
                        project_id=Project.make_gid("Reference"))
        emitter.emit_vertex(p)
    helpers.assert_vertex_file_valid(publication_file)


def test_simple(helpers, emitter_directory, pubmed_file):
    validate(helpers, emitter_directory, pubmed_file)
