import os
import contextlib
import shutil
import pytest

from transform.gtex.cases import transform


@pytest.fixture
def cases_path(request):
    """get the full path of the test input"""
    return os.path.join(request.fspath.dirname, 'source/gtex/GTEx_v7_Annotations_SubjectPhenotypesDS.txt')


@pytest.fixture
def samples_path(request):
    """get the full path of the test input"""
    return os.path.join(request.fspath.dirname, 'source/gtex/GTEx_v7_Annotations_SampleAttributesDS.txt')


def validate(helpers, emitter_directory, cases_path, samples_path):
    """ run xform and test results"""
    aliquot_file = os.path.join(emitter_directory, 'Aliquot.Vertex.json.gz')
    sample_file = os.path.join(emitter_directory, 'Sample.Vertex.json.gz')
    case_file = os.path.join(emitter_directory, 'Case.Vertex.json.gz')
    project_file = os.path.join(emitter_directory, 'Project.Vertex.json.gz')
    program_file = os.path.join(emitter_directory, 'Program.Vertex.json.gz')
    # phenotype_file = os.path.join(emitter_directory, 'Phenotype.Vertex.json.gz')

    programs_edge_file = os.path.join(emitter_directory, 'programs.Edge.json.gz')
    projects_edge_file = os.path.join(emitter_directory, 'projects.Edge.json.gz')
    cases_edge_file = os.path.join(emitter_directory, 'cases.Edge.json.gz')
    case_edge_file = os.path.join(emitter_directory, 'case.Edge.json.gz')
    samples_edge_file = os.path.join(emitter_directory, 'samples.Edge.json.gz')
    sample_edge_file = os.path.join(emitter_directory, 'sample.Edge.json.gz')
    aliquots_edge_file = os.path.join(emitter_directory, 'aliquots.Edge.json.gz')
    # phenotypes_edge_file = os.path.join(emitter_directory, 'phenotypes.Edge.json.gz')

    all_files = [
        # vertices
        aliquot_file, sample_file, case_file, project_file, program_file,
        # edges
        programs_edge_file, projects_edge_file, cases_edge_file, case_edge_file,
        samples_edge_file, sample_edge_file, aliquots_edge_file, aliquots_edge_file
    ]

    # remove output
    with contextlib.suppress(FileNotFoundError):
        shutil.rmtree(emitter_directory)

    # create output
    transform(cases_path=cases_path,
              samples_path=samples_path,
              emitter_prefix=None,
              emitter_directory=emitter_directory)

    for f in all_files:
        if "Vertex.json.gz" in f:
            helpers.assert_vertex_file_valid(f)
        elif "Edge.json.gz" in f:
            helpers.assert_edge_file_valid(f)

    helpers.assert_edge_joins_valid(
        all_files
    )


def test_simple(helpers, emitter_directory, cases_path, samples_path):

    validate(helpers, emitter_directory, cases_path, samples_path)
