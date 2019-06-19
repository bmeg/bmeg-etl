import pytest
import os
import contextlib
import shutil

from transform.pfam.pfam_to_proteins import transform


@pytest.fixture
def data_path(request):
    return os.path.join(request.fspath.dirname, 'source/pfam/homo_sapiens.json')


def test_simple(helpers, emitter_directory, data_path):
    protein_structure_file = os.path.join(emitter_directory, "ProteinStructure.Vertex.json.gz")
    pfam_families_edge_file = os.path.join(emitter_directory, "pfam_families.Edge.json.gz")
    protein_structures_edge_file = os.path.join(emitter_directory, "protein_structures.Edge.json.gz")
    proteins_edge_file = os.path.join(emitter_directory, "proteins.Edge.json.gz")
    protein_edge_file = os.path.join(emitter_directory, "protein.Edge.json.gz")

    all_files = [
        protein_structure_file, pfam_families_edge_file,
        protein_structures_edge_file, proteins_edge_file,
        protein_edge_file
    ]

    # remove output
    with contextlib.suppress(FileNotFoundError):
        shutil.rmtree(emitter_directory)

    # run transform
    transform(data_path=data_path,
              emitter_prefix=None,
              emitter_directory=emitter_directory)

    # check results
    for f in all_files:
        if "Vertex.json.gz" in f:
            helpers.assert_vertex_file_valid(f)
        elif "Edge.json.gz" in f:
            helpers.assert_edge_file_valid(f)
