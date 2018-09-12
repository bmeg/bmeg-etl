
from bmeg.vertex import Compound
from bmeg.emitter import DeduplicationEmitter, new_emitter
import pytest
import os


@pytest.fixture
def emitter_path_prefix(request):
    """ get the full path of the test output """
    return os.path.join(request.fspath.dirname, 'test/test')


def test_dedup_emitter(helpers, emitter_path_prefix):
    """ caching emitter should supress vertexes with same id """
    compounds = [
        Compound(term_id='ANY:foo', term='foo'),
        Compound(term_id='ANY:foo', term='foo'),
        Compound(term_id='ANY:bar', term='bar')
    ]
    emitter = DeduplicationEmitter(prefix=emitter_path_prefix)
    for c in compounds:
        emitter.emit_vertex(c)
    emitter.close()
    compound_file = '{}.Compound.Vertex.json'.format(emitter_path_prefix)
    c = helpers.assert_vertex_file_valid(Compound, compound_file)
    assert c == 2, 'DeduplicationEmitter should have only two vertexes'


def test_new_emitter(emitter_path_prefix):
    emitter = new_emitter('deduplication', prefix=emitter_path_prefix)
    assert emitter


def test_new_emitter_position(emitter_path_prefix):
    emitter = DeduplicationEmitter(prefix=emitter_path_prefix)
    assert emitter


def test_no_edge(emitter_path_prefix):
    """ caching emitter should not accept edge """
    emitter = new_emitter('deduplication', prefix=emitter_path_prefix)
    with pytest.raises(NotImplementedError) as excinfo:
        emitter.emit_edge(None, from_gid=None, to_gid=None)
    assert excinfo, 'DeduplicationEmitter only supports Vertex'
