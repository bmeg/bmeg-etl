""" test maf_transform """

import os
import contextlib
import pytest
import logging
from transform.phenotype.transform import transform
from bmeg.vertex import Phenotype, Aliquot, G2PAssociation
from bmeg.stores import new_store

ALL_FILES = """
normalized.Phenotype.Vertex.json
normalized.PhenotypeOf.Edge.json
normalized.HasPhenotype.Edge.json
""".strip().split()


@pytest.fixture
def output_dir(request):
    """ get the full path of the test output """
    return os.path.join(request.fspath.dirname, 'outputs')


@pytest.fixture
def store_path(request):
    """ get the full path of the test output """
    return os.path.join(request.fspath.dirname, 'source/phenotype/sqlite.db')


def validate(helpers, output_dir, emitter_directory, store_path):

    all_files = [os.path.join(emitter_directory, f) for f in ALL_FILES]
    # remove output
    with contextlib.suppress(FileNotFoundError):
        for f in all_files:
            os.remove(f)
        os.remove(store_path)
    # create output
    transform(output_dir=output_dir, emitter_directory=emitter_directory, store_path=store_path)
    # check output
    compounds = all_files[0]
    phenotype_of = all_files[1]
    has_phenotype = all_files[2]
    helpers.assert_vertex_file_valid(Phenotype, compounds)
    helpers.assert_edge_file_valid(Aliquot, Phenotype, phenotype_of)
    helpers.assert_edge_file_valid(G2PAssociation, Phenotype, has_phenotype)
    # validate vertex for all edges exist
    helpers.assert_edge_joins_valid(all_files, exclude_labels=['Aliquot', 'G2PAssociation'])
    # ensure the store was created
    store = new_store('key-val', path=store_path)
    assert len([c for c in store.all()]) == 19, 'store should have 19 names'


def test_simple(caplog, helpers, output_dir, emitter_directory, store_path):
    """ simple test """
    caplog.set_level(logging.DEBUG)
    validate(helpers, output_dir, emitter_directory, store_path)
