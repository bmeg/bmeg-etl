import pytest
import os
import json
from transform.wiki.compound_transform import transform as compound_transform
from transform.wiki.phenotype_transform import transform as phenotype_transform
from bmeg.ioutils import reader


@pytest.fixture
def phenotype_vertex(request):
    return os.path.join(request.fspath.dirname, 'outputs/phenotype/normalized.Phenotype.Vertex.json')


@pytest.fixture
def output_dir(request):
    return os.path.join(request.fspath.dirname, 'outputs')


@pytest.fixture
def compound_vertex(request):
    return os.path.join(request.fspath.dirname, 'outputs/compound/normalized.Compound.Vertex.json')


def test_phenotype(helpers, output_dir, phenotype_vertex):
    """Simple test to ensure file exists with expected keys."""
    phenotype_transform(vertex_names=phenotype_vertex, output_dir=output_dir, emitter_directory="wiki")
    output_file = f"{output_dir}/wiki/Phenotype.wiki.json.gz"
    assert os.path.exists(output_file)
    with reader(output_file) as ins:
        for line in ins:
            r = json.loads(line)
            assert list(r.keys()) == list(['category', 'title', 'text'])
            if 'Phenotype' == r['category']:
                assert 'SelectedOntology' in r['text']


def test_compound(helpers, output_dir, compound_vertex):
    """Simple test to ensure file exists with expected keys."""
    compound_transform(vertex_names=compound_vertex, output_dir=output_dir, emitter_directory="wiki")
    output_file = f"{output_dir}/wiki/Compound.wiki.json.gz"
    assert os.path.exists(output_file)
    print(output_file)
    with reader(output_file) as ins:
        for line in ins:
            r = json.loads(line)
            assert list(r.keys()) == list(['category', 'title', 'text'])
            if 'Compound' == r['category']:
                assert 'SelectedOntology' in r['text'], r['text']
