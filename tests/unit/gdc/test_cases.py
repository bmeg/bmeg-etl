import os
import shutil
import pytest
import json
from transform.gdc.cases import transform
from bmeg.ioutils import reader


@pytest.fixture
def case_path(request):
    """get the full path of the test input"""
    return os.path.join(request.fspath.dirname, 'source/gdc/cases.json')


def validate(helpers, emitter_directory, case_path):
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
    shutil.rmtree(emitter_directory)

    # create output
    transform(input_path=case_path,
              emitter_directory=emitter_directory)

    for f in all_files:
        if "Vertex.json.gz" in f:
            helpers.assert_vertex_file_valid(f)
        elif "Edge.json.gz" in f:
            helpers.assert_edge_file_valid(f)

    helpers.assert_edge_joins_valid(
        all_files
    )

    # test Aliquot contents
    with reader(aliquot_file) as f:
        for line in f:
            aliquot = json.loads(line)
            assert 'TCGA' not in aliquot['gid'], 'Aliquot gid should be a uuid not {}'.format(aliquot['gid'])
            assert 'TCGA' in aliquot['data']['gdc_attributes']['submitter_id'], 'gdc_attributes.submitter_id should contain TCGA'


def test_simple(helpers, emitter_directory, case_path):

    validate(helpers, emitter_directory, case_path)
