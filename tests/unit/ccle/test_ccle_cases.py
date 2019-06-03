
import os
import contextlib
import pytest
from transform.ccle.ccle_cases import transform
from transform.ccle.depmap_cases import transform as depmap_transform
from bmeg.vertex import Sample, Aliquot, Case, Project, Program


@pytest.fixture
def cellline_meta_path(request):
    return os.path.join(request.fspath.dirname, 'source/ccle/DepMap-2019q1-celllines.csv_v2.csv')


@pytest.fixture
def maf_dir(request):
    return os.path.join(request.fspath.dirname, 'source/ccle/mafs/*')


@pytest.fixture
def expression_path(request):
    return os.path.join(request.fspath.dirname, 'source/ccle/CCLE_depMap_19Q1_TPM.csv')


@pytest.fixture
def drug_response_path(request):
    return os.path.join(request.fspath.dirname, 'source/ccle/CCLE_NP24.2009_Drug_data_2015.02.24.csv')


@pytest.fixture
def cellline_lookup_path(request):
    return os.path.join(request.fspath.dirname, 'source/ccle/cellline_lookup.tsv')


@pytest.fixture
def project_lookup_path(request):
    return os.path.join(request.fspath.dirname, 'source/ccle/cellline_project_lookup.tsv')


@pytest.fixture
def phenotype_lookup_path(request):
    return os.path.join(request.fspath.dirname, 'source/ccle/cellline_phenotype_lookup.tsv')


def validate(helpers, emitter_directory, cellline_meta_path, cellline_lookup_path,
             project_lookup_path, phenotype_lookup_path, drug_response_path,
             expression_path, maf_dir):
    """ run xform and test results"""
    aliquot_file = os.path.join(emitter_directory, 'Aliquot.Vertex.json.gz')
    aliquot_for_file = os.path.join(emitter_directory, 'AliquotFor.Edge.json.gz')
    sample_file = os.path.join(emitter_directory, 'Sample.Vertex.json.gz')
    sample_for_file = os.path.join(emitter_directory, 'SampleFor.Edge.json.gz')
    depmap_case_file = os.path.join(emitter_directory, 'depmap.Case.Vertex.json.gz')
    case_file = os.path.join(emitter_directory, 'Case.Vertex.json.gz')
    in_project_file = os.path.join(emitter_directory, 'InProject.Edge.json.gz')
    project_file = os.path.join(emitter_directory, 'Project.Vertex.json.gz')
    in_program_file = os.path.join(emitter_directory, 'InProgram.Edge.json.gz')
    program_file = os.path.join(emitter_directory, 'Program.Vertex.json.gz')
    phenotype_file = os.path.join(emitter_directory, 'Phenotype.Vertex.json.gz')
    phenotype_of_file = os.path.join(emitter_directory, 'PhenotypeOf.Edge.json.gz')

    all_files = [aliquot_file, aliquot_for_file, sample_file, sample_for_file,
                 depmap_case_file, case_file, in_project_file, project_file,
                 in_program_file, program_file, phenotype_file, phenotype_of_file]

    # remove output
    with contextlib.suppress(FileNotFoundError):
        for f in all_files:
            os.remove(f)

    # create output
    depmap_transform(
        path=cellline_meta_path,
        emitter_prefix="depmap",
        emitter_directory=emitter_directory
    )
    transform(
        cellline_lookup_path=cellline_lookup_path,
        project_lookup_path=project_lookup_path,
        phenotype_lookup_path=phenotype_lookup_path,
        drug_response_path=drug_response_path,
        expression_path=expression_path,
        maf_dir=maf_dir,
        emitter_prefix=None,
        emitter_directory=emitter_directory
    )

    # test.Aliquot.Vertex.json
    helpers.assert_vertex_file_valid(Aliquot, aliquot_file)
    # test.Sample.Vertex.json
    helpers.assert_vertex_file_valid(Sample, sample_file)
    # test.Case.Vertex.json
    depmap_case_count = helpers.assert_vertex_file_valid(Case, depmap_case_file)
    assert depmap_case_count == 10, 'expected case_count'
    # test.Case.Vertex.json
    case_count = helpers.assert_vertex_file_valid(Case, case_file)
    assert case_count == 1, 'expected case_count'
    # test.Project.Vertex.json
    project_count = helpers.assert_vertex_file_valid(Project, project_file)
    assert project_count == 8, 'expected project_count'
    # test.Program.Vertex.json
    program_count = helpers.assert_vertex_file_valid(Program, program_file)
    assert program_count == 1, 'expected program_count'

    # test.AliquotFor.Edge.json
    helpers.assert_edge_file_valid(Aliquot, Sample, aliquot_for_file)

    # test.SampleFor.Edge.json
    helpers.assert_edge_file_valid(Sample, Case, sample_for_file)

    # disabled since InProject is polymorphic
    # test.InProject.Edge.json
    # helpers.assert_edge_file_valid(Case, Project, in_project_file)

    # test.InProgram.Edge.json
    helpers.assert_edge_file_valid(Project, Program, in_program_file)

    # disabled since PhenotypeOf is polymorphic
    # test.PhenotypeOf.Edge.json
    # helpers.assert_edge_file_valid(Case, Phenotype, phenotype_of_file)

    # validate vertex for all edges exist
    helpers.assert_edge_joins_valid(all_files)


def test_simple(helpers, emitter_directory, cellline_meta_path, cellline_lookup_path, project_lookup_path,
                phenotype_lookup_path, drug_response_path, expression_path, maf_dir):

    validate(helpers, emitter_directory, cellline_meta_path, cellline_lookup_path, project_lookup_path,
             phenotype_lookup_path, drug_response_path, expression_path, maf_dir)
