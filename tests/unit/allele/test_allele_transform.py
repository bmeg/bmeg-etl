""" test maf_transform """

import os
import contextlib
import pytest
import logging
import shutil

from transform.allele.harmonize_alleles import transform


@pytest.fixture
def output_directory(request):
    """ get the full path of the test fixture """
    return os.path.join(request.fspath.dirname, 'outputs')


@pytest.fixture
def minimal_allele(request):
    """ get the full path of the test fixture """
    return os.path.join(request.fspath.dirname, 'minimal_alleles.maf')


@pytest.fixture
def annotated_allele(request):
    """ get the full path of the test fixture """
    return os.path.join(request.fspath.dirname, 'annotated_alleles.maf')


def validate(helpers, emitter_directory, output_directory, minimal_allele, annotated_allele):
    allele_file = os.path.join(emitter_directory, 'Allele.Vertex.json.gz')
    allele_gene_file = os.path.join(emitter_directory, 'Allele_Gene_Gene.Edge.json.gz')
    gene_allele_file = os.path.join(emitter_directory, 'Gene_Alleles_Allele.Edge.json.gz')
    allele_pfam_file = os.path.join(emitter_directory, 'Allele_PfamFamily_PfamFamily.Edge.json.gz')
    pfam_allele_file = os.path.join(emitter_directory, 'PfamFamily_Alleles_Allele.Edge.json.gz')
    allele_protein_file = os.path.join(emitter_directory, 'Allele_Protein_Protein.Edge.json.gz')
    protein_allele_file = os.path.join(emitter_directory, 'Protein_Alleles_Allele.Edge.json.gz')
    allele_transcript_file = os.path.join(emitter_directory, 'Allele_Transcript_Transcript.Edge.json.gz')
    transcript_allele_file = os.path.join(emitter_directory, 'Transcript_Alleles_Allele.Edge.json.gz')

    all_files = [allele_file,
                 allele_gene_file, gene_allele_file,
                 allele_pfam_file, pfam_allele_file,
                 allele_protein_file, protein_allele_file,
                 allele_transcript_file, transcript_allele_file]

    # remove output
    with contextlib.suppress(FileNotFoundError):
        shutil.rmtree(emitter_directory)

    # check using memory store
    transform(output_dir=output_directory,
              vertex_filename_pattern='**/*Allele.Vertex.json.gz',
              minimal_maf_file=minimal_allele,
              annotated_maf_file=annotated_allele,
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
        exclude_labels=['Gene', 'Transcript', 'Protein', 'PfamFamily']
    )


def test_simple(caplog, helpers, emitter_directory, output_directory, minimal_allele, annotated_allele):
    """ simple test """
    caplog.set_level(logging.DEBUG)
    validate(helpers, emitter_directory, output_directory, minimal_allele, annotated_allele)
