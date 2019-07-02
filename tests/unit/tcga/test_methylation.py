import os
import contextlib
import pytest
import shutil
from transform.tcga.methylation import transform


@pytest.fixture
def source_path(request):
    """ get the full path of the test input """
    return os.path.join(request.fspath.dirname, 'source/tcga/methylation/IlluminaHumanMethylation450.tar.gz')


@pytest.fixture
def id_lookup_path(request):
    """ get the full path of the test fixture """
    return os.path.join(request.fspath.dirname, 'source/gdc/id_lookup.tsv')


@pytest.fixture
def project_lookup_path(request):
    """ get the full path of the test fixture """
    return os.path.join(request.fspath.dirname, 'source/gdc/project_lookup.tsv')


def validate(helpers, emitter_directory, source_path, id_lookup_path, project_lookup_path):
    """ run xform and test results"""
    methylation_probe_file = os.path.join(emitter_directory, 'IlluminaHumanMethylation450.MethylationProbe.Vertex.json.gz')
    methylation_file = os.path.join(emitter_directory, 'IlluminaHumanMethylation450.Methylation.Vertex.json.gz')

    methylation_probes_edge_file = os.path.join(emitter_directory, 'IlluminaHumanMethylation450.MethylationProbe_Gene_Gene.Edge.json.gz')
    genes_edge_file = os.path.join(emitter_directory, 'IlluminaHumanMethylation450.Gene_MethylationProbes_MethylationProbe.Edge.json.gz')
    methylations_edge_file = os.path.join(emitter_directory, 'IlluminaHumanMethylation450.Methylation_Aliquot_Aliquot.Edge.json.gz')
    aliquot_edge_file = os.path.join(emitter_directory, 'IlluminaHumanMethylation450.Aliquot_Methylations_Methylation.Edge.json.gz')

    all_files = [methylation_probe_file, methylation_file,
                 methylation_probes_edge_file, genes_edge_file,
                 methylations_edge_file, aliquot_edge_file]

    # remove output
    with contextlib.suppress(FileNotFoundError):
        shutil.rmtree(emitter_directory)

    # create output
    transform(source_path=source_path,
              id_lookup_path=id_lookup_path,
              project_lookup_path=project_lookup_path,
              emitter_prefix=None,
              emitter_directory=emitter_directory)

    # ratify
    for f in all_files:
        if "Vertex.json.gz" in f:
            helpers.assert_vertex_file_valid(f)
        elif "Edge.json.gz" in f:
            helpers.assert_edge_file_valid(f)

    helpers.assert_edge_joins_valid(
        [methylation_file, methylations_edge_file, aliquot_edge_file],
        exclude_labels=['Aliquot']
    )


def test_simple(helpers, emitter_directory, source_path, id_lookup_path, project_lookup_path):
    """ just run validate"""
    validate(helpers, emitter_directory, source_path, id_lookup_path, project_lookup_path)
