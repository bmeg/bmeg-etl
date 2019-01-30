import pytest
import os
from transform.ensembl.proteins import transform
from bmeg.vertex import Protein


@pytest.fixture
def protein_table_path(request):
    """ get the full path of the test output """
    return os.path.join(request.fspath.dirname, 'source/ensembl/Homo_sapiens.GRCh37.85.uniprot.tsv.gz')


def test_simple(helpers, emitter_directory, protein_table_path):
    protein_file = os.path.join(emitter_directory, 'Protein.Vertex.json.gz')
    proteinfor_file = os.path.join(emitter_directory, 'ProteinFor.Edge.json.gz')
    transform(protein_table_path=protein_table_path, emitter_directory=emitter_directory)
    helpers.assert_vertex_file_valid(Protein, '{}/Protein.Vertex.json.gz'.format(emitter_directory))
    helpers.assert_edge_joins_valid([protein_file, proteinfor_file], exclude_labels=['Transcript'])
