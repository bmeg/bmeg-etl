import pytest
import os
import contextlib
import shutil

from transform.ensembl.uniprot import transform


@pytest.fixture
def protein_table_path(request):
    """ get the full path of the test output """
    return os.path.join(request.fspath.dirname, 'source/ensembl/Homo_sapiens.GRCh37.85.uniprot.tsv.gz')


def test_simple(helpers, emitter_directory, protein_table_path):
    uniprot_file = os.path.join(emitter_directory, 'Uniprot.Vertex.json.gz')
    tpp_edge_file = os.path.join(emitter_directory, 'Protein_Uniprot_Uniprot.Edge.json.gz')
    ptt_edge_file = os.path.join(emitter_directory, 'Uniprot_Protein_Protein.Edge.json.gz')

    # remove output
    with contextlib.suppress(FileNotFoundError):
        shutil.rmtree(emitter_directory)

    transform(protein_table_path=protein_table_path, emitter_directory=emitter_directory)

    helpers.assert_vertex_file_valid(uniprot_file)
    helpers.assert_edge_file_valid(tpp_edge_file)
    helpers.assert_edge_file_valid(ptt_edge_file)
    helpers.assert_edge_joins_valid([uniprot_file, tpp_edge_file, ptt_edge_file],
                                    exclude_labels=['Protein'])
