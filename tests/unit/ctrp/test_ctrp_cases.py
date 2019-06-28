import contextlib
import os
import shutil
import pytest
from transform.ctrp.cases import transform
from transform.ccle.depmap_cases import transform as depmap_transform


@pytest.fixture
def metadrugPath(request):
    return os.path.join(request.fspath.dirname, 'source/ctrp/v20.meta.per_compound.txt')


@pytest.fixture
def metacelllinePath(request):
    return os.path.join(request.fspath.dirname, 'source/ctrp/v20.meta.per_cell_line.txt')


@pytest.fixture
def responsePath(request):
    return os.path.join(request.fspath.dirname, 'source/ctrp/v20.data.curves_post_qc.txt')


@pytest.fixture
def metaexperimentPath(request):
    return os.path.join(request.fspath.dirname, 'source/ctrp/v20.meta.per_experiment.txt')


@pytest.fixture
def curvePath(request):
    return os.path.join(request.fspath.dirname, 'source/ctrp/v20.data.per_cpd_post_qc.txt')


@pytest.fixture
def cellline_meta_path(request):
    return os.path.join(request.fspath.dirname, 'source/ccle/DepMap-2019q1-celllines.csv_v2.csv')


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
             project_lookup_path, phenotype_lookup_path, metadrugPath,
             metacelllinePath, responsePath, metaexperimentPath, curvePath):
    """ run xform and test results"""
    aliquot_file = os.path.join(emitter_directory, 'Aliquot.Vertex.json.gz')
    sample_file = os.path.join(emitter_directory, 'Sample.Vertex.json.gz')
    case_file = os.path.join(emitter_directory, 'Case.Vertex.json.gz')
    project_file = os.path.join(emitter_directory, 'Project.Vertex.json.gz')
    program_file = os.path.join(emitter_directory, 'Program.Vertex.json.gz')
    phenotype_file = os.path.join(emitter_directory, 'Phenotype.Vertex.json.gz')

    proj_pp_edge_file = os.path.join(emitter_directory, 'Project_Programs_Program.Edge.json.gz')
    prog_pps_edge_file = os.path.join(emitter_directory, 'Program_Projects_Project.Edge.json.gz')
    cpp_edge_file = os.path.join(emitter_directory, 'Case_Projects_Project.Edge.json.gz')
    pcc_edge_file = os.path.join(emitter_directory, 'Project_Cases_Case.Edge.json.gz')
    scc_edge_file = os.path.join(emitter_directory, 'Sample_Case_Case.Edge.json.gz')
    css_edge_file = os.path.join(emitter_directory, 'Case_Samples_Sample.Edge.json.gz')
    ass_edge_file = os.path.join(emitter_directory, 'Aliquot_Samples_Sample.Edge.json.gz')
    saa_edge_file = os.path.join(emitter_directory, 'Sample_Aliquots_Aliquot.Edge.json.gz')
    case_pp_edge_file = os.path.join(emitter_directory, 'Case_Phenotypes_Phenotype.Edge.json.gz')
    pheno_cc_edge_file = os.path.join(emitter_directory, 'Phenotype_Cases_Case.Edge.json.gz')
    sample_pp_edge_file = os.path.join(emitter_directory, 'Sample_Phenotypes_Phenotype.Edge.json.gz')
    pheno_sample_edge_file = os.path.join(emitter_directory, 'Phenotype_Samples_Sample.Edge.json.gz')

    all_files = [
        # vertices
        aliquot_file, sample_file, case_file, project_file,
        program_file, phenotype_file,
        # edges
        proj_pp_edge_file, prog_pps_edge_file, cpp_edge_file, pcc_edge_file,
        scc_edge_file, css_edge_file, ass_edge_file, saa_edge_file, case_pp_edge_file,
        pheno_cc_edge_file, sample_pp_edge_file, pheno_sample_edge_file
    ]

    # remove output
    with contextlib.suppress(FileNotFoundError):
        shutil.rmtree(emitter_directory)

    # create output
    depmap_transform(
        path=cellline_meta_path,
        phenotype_lookup_path=phenotype_lookup_path,
        emitter_prefix="depmap",
        emitter_directory=emitter_directory
    )
    transform(
        cellline_lookup_path=cellline_lookup_path,
        project_lookup_path=project_lookup_path,
        phenotype_lookup_path=phenotype_lookup_path,
        metadrugPath=metadrugPath,
        metacelllinePath=metacelllinePath,
        responsePath=responsePath,
        metaexperimentPath=metaexperimentPath,
        curvePath=curvePath,
        emitter_prefix=None,
        emitter_directory=emitter_directory
    )

    for f in all_files:
        if "Vertex.json.gz" in f:
            helpers.assert_vertex_file_valid(f)
        elif "Edge.json.gz" in f:
            helpers.assert_edge_file_valid(f)

    helpers.assert_edge_joins_valid(
        all_files,
        exclude_labels=["DrugResponse"]
    )


def test_simple(helpers, emitter_directory, cellline_meta_path, cellline_lookup_path,
                project_lookup_path, phenotype_lookup_path, metadrugPath,
                metacelllinePath, responsePath, metaexperimentPath, curvePath):

    validate(helpers, emitter_directory, cellline_meta_path, cellline_lookup_path,
             project_lookup_path, phenotype_lookup_path, metadrugPath,
             metacelllinePath, responsePath, metaexperimentPath, curvePath)
