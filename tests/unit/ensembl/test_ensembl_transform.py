import pytest
import os
import contextlib
import shutil

from transform.ensembl.transform import transform


@pytest.fixture
def gff3_path(request):
    """ get the full path of the test output """
    return os.path.join(request.fspath.dirname, 'source/ensembl/Homo_sapiens.GRCh37.87.gff3.gz')


def test_simple(emitter_directory, gff3_path, helpers):
    exon_file = '{}/Exon.Vertex.json.gz'.format(emitter_directory)
    transcript_file = '{}/Transcript.Vertex.json.gz'.format(emitter_directory)
    gene_file = '{}/Gene.Vertex.json.gz'.format(emitter_directory)
    protein_file = '{}/Protein.Vertex.json.gz'.format(emitter_directory)
    tee_edge_file = '{}/Transcript_Exons_Exon.Edge.json.gz'.format(emitter_directory)
    ett_edge_file = '{}/Exon_Transcripts_Transcript.Edge.json.gz'.format(emitter_directory)
    gtt_edge_file = '{}/Gene_Transcripts_Transcript.Edge.json.gz'.format(emitter_directory)
    tgg_edge_file = '{}/Transcript_Gene_Gene.Edge.json.gz'.format(emitter_directory)
    tpp_edge_file = os.path.join(emitter_directory, 'Transcript_Protein_Protein.Edge.json.gz')
    ptt_edge_file = os.path.join(emitter_directory, 'Protein_Transcript_Transcript.Edge.json.gz')

    # remove output
    with contextlib.suppress(FileNotFoundError):
        shutil.rmtree(emitter_directory)

    transform(gff3_path=gff3_path, emitter_directory=emitter_directory)

    exon_count = helpers.assert_vertex_file_valid(exon_file)
    transcript_count = helpers.assert_vertex_file_valid(transcript_file)
    gene_count = helpers.assert_vertex_file_valid(gene_file)
    protein_count = helpers.assert_vertex_file_valid(protein_file)
    assert exon_count == 624  # there are 843 entry lines, but 624 unique ids
    assert transcript_count == 211
    assert gene_count == 81
    assert protein_count == 38
    helpers.assert_edge_file_valid(tee_edge_file)
    helpers.assert_edge_file_valid(ett_edge_file)
    helpers.assert_edge_file_valid(tgg_edge_file)
    helpers.assert_edge_file_valid(gtt_edge_file)
    helpers.assert_edge_file_valid(tpp_edge_file)
    helpers.assert_edge_file_valid(ptt_edge_file)
    helpers.assert_edge_joins_valid([
        exon_file, transcript_file, gene_file, protein_file,
        tee_edge_file, ett_edge_file, gtt_edge_file, tgg_edge_file,
        tpp_edge_file, ptt_edge_file
    ])
