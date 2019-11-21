""" test maf_transform """

import os
import contextlib
import pytest
import logging
import json
import shutil

from transform.allele.harmonize_alleles import create_minimal_maf, run_maf2maf, transform
from bmeg.ioutils import reader


@pytest.fixture
def output_directory(request):
    """ get the full path of the test fixture """
    return os.path.join(request.fspath.dirname, 'outputs')


def validate(helpers, emitter_directory, output_directory):
    allele_file = os.path.join(emitter_directory, 'Allele.Vertex.json.gz')
    allele_gene_file = os.path.join(emitter_directory, 'Allele_Gene_Gene.Edge.json.gz')
    allele_pfam_file = os.path.join(emitter_directory, 'Allele_PfamFamily_PfamFamily.Edge.json.gz')
    allele_protein_file = os.path.join(emitter_directory, 'Allele_Protein_Protein.Edge.json.gz')
    allele_transcript_file = os.path.join(emitter_directory, 'Allele_Transcript_Transcript.Edge.json.gz')
    allele_somaticcallset_file = os.path.join(emitter_directory, 'Allele_SomaticCallset_SomaticCallset.Edge.json.gz')
    # backrefs TODO

    all_files = [allele_file, allele_gene_file, allele_pfam_file, allele_protein_file,
                 allele_transcript_file, allele_somaticcallset_file]

    # remove output
    with contextlib.suppress(FileNotFoundError):
        shutil.rmtree(emitter_directory)

    # check using memory store
    transform(output_dir=output_directory,
              vertex_filename_pattern='**/*Allele.Vertex.json.gz',
              minimal_maf_file=os.path.join(emitter_directory, 'minimal_alleles.maf'),
              annotated_maf_file=os.path.join(emitter_directory, 'annotated_alleles.maf'),
              emitter_directory=emitter_directory,
              emitter_prefix=None)

    for f in all_files:
        if "Vertex.json.gz" in f:
            count = helpers.assert_vertex_file_valid(f)
            if f == allele_file:
                assert count == 30
        elif "Edge.json.gz" in f:
            helpers.assert_edge_file_valid(f)

    # validate vertex for all edges exist
    helpers.assert_edge_joins_valid(
        all_files,
        exclude_labels=['SomaticCallset']
    )


def test_simple(caplog, helpers, emitter_directory, output_directory):
    """ simple test """
    caplog.set_level(logging.DEBUG)
    validate(helpers, emitter_directory, output_directory)
