import os
import contextlib
import pytest
import shutil

from transform.publication.transform import transform


@pytest.fixture
def output_dir(request):
    """ get the full path of the test fixture """
    return os.path.join(request.fspath.dirname, 'outputs')


def validate(helpers, emitter_directory, output_dir):
    publication_file = os.path.join(emitter_directory, 'stub.Publication.Vertex.json.gz')
    publications_edge_file = os.path.join(output_dir, 'g2p/G2PAssociation_Publications_Publication.Edge.json.gz')

    # remove output
    with contextlib.suppress(FileNotFoundError):
        shutil.rmtree(emitter_directory)

    transform(
        output_dir=output_dir,
        emitter_directory=emitter_directory
    )

    helpers.assert_vertex_file_valid(publication_file)
    helpers.assert_edge_joins_valid(
        [publication_file, publications_edge_file],
        exclude_labels=["G2PAssociation"]
    )


def test_simple(helpers, emitter_directory, output_dir):
    validate(helpers, emitter_directory, output_dir)
