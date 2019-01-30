import pytest
import os
from transform.ensembl.transform import transform
from bmeg.vertex import Transcript, Exon, Gene


@pytest.fixture
def gff3_path(request):
    """ get the full path of the test output """
    return os.path.join(request.fspath.dirname, 'source/ensembl/Homo_sapiens.GRCh37.87.gff3.gz')


def test_simple(emitter_directory, gff3_path, helpers):
    exon_file = '{}/Exon.Vertex.json.gz'.format(emitter_directory)
    exonfor_file = '{}/ExonFor.Edge.json.gz'.format(emitter_directory)
    transcript_file = '{}/Transcript.Vertex.json.gz'.format(emitter_directory)
    transcriptfor_file = '{}/TranscriptFor.Edge.json.gz'.format(emitter_directory)
    gene_file = '{}/Gene.Vertex.json.gz'.format(emitter_directory)

    transform(gff3_path=gff3_path, emitter_directory=emitter_directory)

    exon_count = helpers.assert_vertex_file_valid(Exon, exon_file)
    transcript_count = helpers.assert_vertex_file_valid(Transcript, transcript_file)
    gene_count = helpers.assert_vertex_file_valid(Gene, gene_file)
    assert exon_count == 624  # there are 843 entry lines, but 624 unique ids
    assert transcript_count == 211
    assert gene_count == 81
    helpers.assert_edge_joins_valid([exon_file, transcript_file, gene_file, exonfor_file, transcriptfor_file])
