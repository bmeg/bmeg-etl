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
    protein_pfam_edge_file = os.path.join(emitter_directory, "Protein_PfamFamilies_PfamFamily.Edge.json.gz")
    pfam_protein_edge_file = os.path.join(emitter_directory, "PfamFamily_Proteins_Protein.Edge.json.gz")
    protein_struct_edge_file = os.path.join(emitter_directory, "Protein_ProteinStructures_ProteinStructure.Edge.json.gz")
    struct_protein_edge_file = os.path.join(emitter_directory, "ProteinStructure_Protein_Protein.Edge.json.gz")

    all_files = [
        protein_structure_file, protein_pfam_edge_file, pfam_protein_edge_file,
        protein_struct_edge_file, struct_protein_edge_file
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
