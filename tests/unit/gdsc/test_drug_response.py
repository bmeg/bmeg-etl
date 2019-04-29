
import os
import contextlib
import pytest
from transform.gdsc.drug_response import transform
from bmeg.vertex import Compound, DrugResponse, Aliquot, Project


@pytest.fixture
def GDSC_AUC_file(request):
    """ get the full path of the test fixture """
    return os.path.join(request.fspath.dirname, 'source/gdsc/GDSC_AUC.csv')


@pytest.fixture
def GDSC_IC50_file(request):
    """ get the full path of the test fixture """
    return os.path.join(request.fspath.dirname, 'source/gdsc/GDSC_IC50.csv')


@pytest.fixture
def metadrugPath(request):
    return os.path.join(request.fspath.dirname, 'source/gdsc/Screened_Compounds.xlsx')


@pytest.fixture
def project_lookup_path(request):
    return os.path.join(request.fspath.dirname, 'source/ccle/cellline_project_lookup.tsv')


def validate(helpers, emitter_directory, metadrugPath, GDSC_AUC_file, GDSC_IC50_file, project_lookup_path):

    profile_file = os.path.join(emitter_directory, 'DrugResponse.Vertex.json.gz')
    profile_in_file = os.path.join(emitter_directory, 'ResponseIn.Edge.json.gz')
    response_to_file = os.path.join(emitter_directory, 'ResponseTo.Edge.json.gz')
    compound_file = os.path.join(emitter_directory, 'Compound.Vertex.json.gz')
    tested_in_file = os.path.join(emitter_directory, 'TestedIn.Edge.json.gz')

    all_files = [profile_file, profile_in_file, response_to_file, compound_file, tested_in_file]

    # remove output
    with contextlib.suppress(FileNotFoundError):
        for f in all_files:
            os.remove(f)

    # create output
    transform(
        project_lookup_path=project_lookup_path,
        drugs_meta_path=metadrugPath,
        auc_path=GDSC_AUC_file,
        ic50_path=GDSC_IC50_file,
        emitter_prefix=None,
        emitter_directory=emitter_directory
    )

    # ratify
    helpers.assert_vertex_file_valid(DrugResponse, profile_file)
    helpers.assert_edge_file_valid(DrugResponse, Aliquot, profile_in_file)
    helpers.assert_edge_file_valid(DrugResponse, Compound, response_to_file)
    helpers.assert_vertex_file_valid(Compound, compound_file)
    helpers.assert_edge_file_valid(Compound, Project, tested_in_file)
    helpers.assert_edge_joins_valid(all_files, exclude_labels=['Aliquot', 'Project'])


def test_simple(helpers, emitter_directory, metadrugPath, GDSC_AUC_file, GDSC_IC50_file, project_lookup_path):
    validate(helpers, emitter_directory, metadrugPath, GDSC_AUC_file, GDSC_IC50_file, project_lookup_path)
