import pytest
import os
from transform.ensembl.transform import transform
from bmeg.vertex import Transcript, Exon, Gene


@pytest.fixture
def gff3_path(request):
    """ get the full path of the test output """
    return os.path.join(request.fspath.dirname, 'source/ensembl/Homo_sapiens.GRCh37.87.gff3.gz')


def test_simple(emitter_directory, gff3_path, helpers):
    """ get the missing transcripts """
    transform(gff3_path=gff3_path, emitter_directory=emitter_directory)
    transcript_count = helpers.assert_vertex_file_valid(Transcript, '{}/Transcript.Vertex.json.gz'.format(emitter_directory))
    exon_count = helpers.assert_vertex_file_valid(Exon, '{}/Exon.Vertex.json.gz'.format(emitter_directory))
    gene_count = helpers.assert_vertex_file_valid(Gene, '{}/Gene.Vertex.json.gz'.format(emitter_directory))
    assert transcript_count == 211
    assert exon_count == 624 # there are 843 entry lines, but 624 unique ids
    assert gene_count == 81
