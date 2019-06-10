import pytest
import os
from transform.ensembl.proteins import transform


@pytest.fixture
def protein_table_path(request):
    """ get the full path of the test output """
    return os.path.join(request.fspath.dirname, 'source/ensembl/Homo_sapiens.GRCh37.85.uniprot.tsv.gz')


def test_simple(helpers, emitter_directory, protein_table_path):
    protein_file = os.path.join(emitter_directory, 'Protein.Vertex.json.gz')
    transcript_edge_file = os.path.join(emitter_directory, 'transcript.Edge.json.gz')
    protein_edge_file = os.path.join(emitter_directory, 'protein.Edge.json.gz')

    transform(protein_table_path=protein_table_path, emitter_directory=emitter_directory)

    helpers.assert_vertex_file_valid(protein_file)
    helpers.assert_edge_file_valid(transcript_edge_file)
    helpers.assert_edge_file_valid(protein_edge_file)
    helpers.assert_edge_joins_valid([protein_file, protein_edge_file, transcript_edge_file], exclude_labels=['Transcript'])
