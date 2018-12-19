import os
import contextlib
import pytest
from transform.tcga.methylation import transform
from bmeg.vertex import Methylation, MethylationProbe, Aliquot, Gene


@pytest.fixture
def source_path(request):
    """ get the full path of the test input """
    return os.path.join(request.fspath.dirname, 'source/tcga/methylation/IlluminaHumanMethylation450.tar.gz')


@pytest.fixture
def aliquot_source_path(request):
    """ get the full path of the test aliquot input """
    return os.path.join(request.fspath.dirname, 'source/tcga/methylation/Aliquot.Vertex.json.gz')


def validate(helpers, source_path, aliquot_source_path, emitter_directory):
    """ run xform and test results"""
    kvstore_path = os.path.join(emitter_directory, 'methylation.db')
    methylation_probe_file = os.path.join(emitter_directory, 'IlluminaHumanMethylation450.MethylationProbe.Vertex.json.gz')
    methylation_probe_for_file = os.path.join(emitter_directory, 'IlluminaHumanMethylation450.MethylationProbeFor.Edge.json.gz')
    methylation_file = os.path.join(emitter_directory, 'IlluminaHumanMethylation450.Methylation.Vertex.json.gz')
    methylation_of_file = os.path.join(emitter_directory, 'IlluminaHumanMethylation450.MethylationOf.Edge.json.gz')

    all_files = [methylation_probe_file, methylation_probe_for_file, methylation_file, methylation_of_file, kvstore_path]
    # remove output
    with contextlib.suppress(FileNotFoundError):
        for f in all_files:
            os.remove(f)

    # create output
    transform(source_path=source_path,
              aliquot_path=aliquot_source_path,
              kvstore_path=kvstore_path,
              emitter_directory=emitter_directory)

    # ratify
    helpers.assert_vertex_file_valid(Methylation, methylation_file)
    helpers.assert_edge_file_valid(Methylation, Aliquot, methylation_of_file)
    helpers.assert_vertex_file_valid(MethylationProbe, methylation_probe_file)
    helpers.assert_edge_file_valid(MethylationProbe, Gene, methylation_probe_for_file)
    all_files = [methylation_probe_file, methylation_probe_for_file, methylation_file, methylation_of_file]
    helpers.assert_edge_joins_valid(all_files, exclude_labels=['Aliquot', 'Gene', 'MethylationProbe'])


def test_simple(helpers, source_path, aliquot_source_path, emitter_directory):
    """ just run validate"""
    validate(helpers, source_path, aliquot_source_path, emitter_directory)
