""" test maf_transform """

import os
import contextlib
import pytest
import logging
import shutil
from transform.phenotype.transform import transform
from bmeg.stores import new_store


ALL_FILES = """
normalized.Phenotype.Vertex.json.gz
normalized.Sample_Phenotypes_Phenotype.Edge.json.gz
normalized.G2PAssociation_Phenotypes_Phenotype.Edge.json.gz
""".strip().split()


@pytest.fixture
def output_dir(request):
    """ get the full path of the test output """
    return os.path.join(request.fspath.dirname, 'outputs')


@pytest.fixture
def store_path(request):
    """ get the full path of the test output """
    return os.path.join(request.fspath.dirname, 'source/phenotype/sqlite.db')


def validate(helpers, emitter_directory, output_dir, store_path):

    all_files = [os.path.join(emitter_directory, f) for f in ALL_FILES]

    # remove output
    with contextlib.suppress(FileNotFoundError):
        shutil.rmtree(emitter_directory)
        os.remove(store_path)

    # create output
    transform(output_dir=output_dir,
              emitter_directory=emitter_directory,
              store_path=store_path)

    # check output
    phenotypes = all_files[0]
    phenotype_of = all_files[1]
    has_phenotype = all_files[2]
    helpers.assert_vertex_file_valid(phenotypes)
    helpers.assert_edge_file_valid(phenotype_of)
    helpers.assert_edge_file_valid(has_phenotype)
    # validate vertex for all edges exist
    helpers.assert_edge_joins_valid(all_files, exclude_labels=['Sample', 'G2PAssociation'])

    # ensure the store was created
    store = new_store('key-val', path=store_path)
    assert len([c for c in store.all()]) == 19, 'store should have 19 names'


def test_simple(caplog, helpers, emitter_directory, output_dir, store_path):
    """ simple test """
    caplog.set_level(logging.DEBUG)
    validate(helpers, emitter_directory, output_dir, store_path)
