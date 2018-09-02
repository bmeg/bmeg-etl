""" test maf_transform """

import os
import contextlib
import pytest
import logging
import json
from transform.allele.transform import transform, sort_allele_files, group_sorted_alleles, merge
from bmeg.vertex import Allele, AlleleAnnotations
from bmeg.ioutils import reader

@pytest.fixture
def output_directory(request):
    """ get the full path of the test fixture """
    return os.path.join(request.fspath.dirname, 'outputs')


@pytest.fixture
def emitter_path_prefix(request):
    """ get the full path of the test output """
    return os.path.join(request.fspath.dirname, 'test/test')


@pytest.fixture
def myvariantinfo_path(request):
    """ get the full path of the test input """
    return os.path.join(request.fspath.dirname, 'source/myvariantinfo/biothings_current_old_hg19.json.gz')


def validate_myvariantinfo_count(allele_file):
    myvariantinfo_count = 0
    with open(allele_file, 'r', encoding='utf-8') as f:
        for line in f:
            # should be json
            allele = json.loads(line)
            # optional keys, if set should be non null
            if allele['data']['annotations']['myvariantinfo']:
                myvariantinfo_count += 1
    assert myvariantinfo_count == 45, 'we should have had myvariantinfo hits'


def validate(helpers, output_directory, emitter_path_prefix, myvariantinfo_path):
    allele_file = '{}.Allele.Vertex.json'.format(emitter_path_prefix)
    # remove output
    with contextlib.suppress(FileNotFoundError):
        os.remove(allele_file)

    # check using memory store
    transform(output_directory,
              prefix=emitter_path_prefix,
              vertex_filename_pattern='**/*.Allele.Vertex.json',
              myvariant_store_name='myvariantinfo-memory',
              myvariantinfo_path=myvariantinfo_path)
    # test/test.Allele.Vertex.json
    helpers.assert_vertex_file_valid(Allele, allele_file)
    validate_myvariantinfo_count(allele_file)

    # remove output
    with contextlib.suppress(FileNotFoundError):
        os.remove(allele_file)

    # check using sqllite
    transform(output_directory,
              prefix=emitter_path_prefix,
              vertex_filename_pattern='**/*.Allele.Vertex.json',
              myvariant_store_name='myvariantinfo-sqlite',
              myvariantinfo_path=myvariantinfo_path)
    # test/test.Allele.Vertex.json
    helpers.assert_vertex_file_valid(Allele, allele_file)
    validate_myvariantinfo_count(allele_file)


def test_simple(caplog, helpers, output_directory, emitter_path_prefix, myvariantinfo_path):
    """ simple test """
    caplog.set_level(logging.DEBUG)
    validate(helpers, output_directory, emitter_path_prefix, myvariantinfo_path)


def test_sort_allele_files(output_directory):
    """ ensure that allele files are sorted """
    path = '{}/{}'.format(output_directory, '**/*.Allele.Vertex.json')
    sorted_allele_file = sort_allele_files(path, '/tmp')
    with reader(sorted_allele_file) as ins:
        _id = 'zzzzzzzzzzz'
        for line in ins:
            data = json.loads(line)
            assert data['_id'] < _id or data['_id'] == _id
            _id = data['_id']


def test_group_sorted_alleles(output_directory):
    """ ensure that allele files are sorted """
    path = '{}/{}'.format(output_directory, '**/*.Allele.Vertex.json')
    sorted_allele_file = sort_allele_files(path, '/tmp')
    t = 0
    uniq_ids = []
    for alleles in group_sorted_alleles(sorted_allele_file):
        assert len(alleles) > 0, 'should be at least one allele'
        ids = set([allele.gid() for allele in alleles])
        assert len(ids) == 1, 'all alleles should have the same id'
        t += len(alleles)
        _id = list(ids)[0]
        assert _id not in uniq_ids, 'we should not have seen this id before'
        uniq_ids.append(_id)
    assert t == 274, 'there should be 274 alleles'


def test_merge():
    """ merge works """
    alleles = [
        Allele(genome='G', chromosome='1', start=1, end=2, reference_bases='A', alternate_bases='G', annotations=AlleleAnnotations()),
        Allele(genome='G', chromosome='1', start=1, end=2, reference_bases='A', alternate_bases='G', annotations=AlleleAnnotations(maf={'a': "A"})),
        Allele(genome='G', chromosome='1', start=1, end=2, reference_bases='A', alternate_bases='G', annotations=AlleleAnnotations(mc3={'a': "B"})),
        Allele(genome='G', chromosome='1', start=1, end=2, reference_bases='A', alternate_bases='G', annotations=AlleleAnnotations(ccle={'a': "C"})),
        Allele(genome='G', chromosome='1', start=1, end=2, reference_bases='A', alternate_bases='G', annotations=AlleleAnnotations(myvariantinfo={'a': "D"})),
    ]
    allele = merge(alleles)
    assert allele, 'should return an allele'
    assert allele.annotations, 'should return an allele.annotations'
    assert allele.annotations.maf, 'should return an allele.annotations.maf'
    assert allele.annotations.mc3, 'should return an allele.annotations.mc3'
    assert allele.annotations.ccle, 'should return an allele.annotations.ccle'
    assert allele.annotations.myvariantinfo, 'should return an allele.annotations.myvariantinfo'
