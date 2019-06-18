import pytest
import os
import shutil

from transform.ensembl.missing_transcripts import transform


@pytest.fixture
def missing_transcript_ids_filename(request):
    """ get the full path of the test output """
    return os.path.join(request.fspath.dirname, 'source/ensembl/missing_transcript_ids.txt')


def test_simple(emitter_directory, missing_transcript_ids_filename, helpers):
    """ get the missing transcripts """

    # remove output
    shutil.rmtree(emitter_directory)

    transform(output_dir=emitter_directory,
              prefix='missing',
              missing_transcript_ids_filename=missing_transcript_ids_filename)

    transcript_count = helpers.assert_vertex_file_valid('{}/missing.Transcript.Vertex.json.gz'.format(emitter_directory))
    assert transcript_count == 3
