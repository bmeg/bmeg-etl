
from bmeg.vertex import Compound
from bmeg.emitter import DeduplicationEmitter, new_emitter
import pytest


def test_dedup_emitter(helpers, emitter_directory, emitter_prefix):
    """ caching emitter should supress vertexes with same id """
    compounds = [
        Compound(term_id='ANY:foo', term='foo'),
        Compound(term_id='ANY:foo', term='foo'),
        Compound(term_id='ANY:bar', term='bar')
    ]
    emitter = DeduplicationEmitter(prefix=emitter_prefix, directory=emitter_directory)
    for c in compounds:
        emitter.emit_vertex(c)
    emitter.close()
    compound_file = '{}/{}.Compound.Vertex.json'.format(emitter_directory, emitter_prefix)
    c = helpers.assert_vertex_file_valid(Compound, compound_file)
    assert c == 2, 'DeduplicationEmitter should have only two vertexes'


def test_new_emitter(emitter_directory, emitter_prefix):
    emitter = new_emitter('deduplication', prefix=emitter_prefix, directory=emitter_directory)
    assert emitter


def test_no_edge(emitter_directory, emitter_prefix):
    """ caching emitter should not accept edge """
    emitter = new_emitter('deduplication', prefix=emitter_prefix, directory=emitter_directory)
    with pytest.raises(NotImplementedError) as excinfo:
        emitter.emit_edge(None, from_gid=None, to_gid=None)
    assert excinfo, 'DeduplicationEmitter only supports Vertex'
