import contextlib
from bmeg import Compound
from bmeg.stores import new_store
from bmeg.enrichers.drug_enricher import compound_factory


def validate_put_get(s):
    s.put('a', 'AAAA')
    a = s.get('a')
    assert a == 'AAAA', 'should get what we put'


def validate_batch(s):
    batch = [('a', 'AAAA'), ('b', 'BBBB')]
    s.load_many(batch)
    batch2 = [b for b in s.all()]
    assert ['AAAA', 'BBBB'] == batch2, 'should get all what we loaded'
    ids = [k for k in s.all_ids()]
    assert ['a', 'b'] == ids, 'should get all ids loaded'


def validate(s):
    validate_put_get(s)
    validate_batch(s)
    assert s.backend(), 'should return implementation backend'


def validate_dataclass(s):
    c = compound_factory('foo')
    c_id = c.gid()
    s.put(c)
    c = s.get(c_id)
    assert c.submitter_id == 'foo', 'should get what we put'
    batch = [compound_factory('foo'), compound_factory('bar')]
    s.load_many(batch)
    batch2 = [b.submitter_id for b in s.all()]
    assert ['foo', 'bar'] == batch2, 'should get all what we loaded'
    ids = sorted([k for k in s.all_ids()])
    assert sorted(['Compound:TODO:foo', 'Compound:TODO:bar']) == ids, 'should get all ids loaded'
    assert s.backend(), 'should return implementation backend'


def test_memory():
    validate(new_store('memory'))


def test_keyval():
    validate(new_store('key-val', index=True))


def test_keyval_no_index():
    validate_put_get(new_store('key-val', index=False))


def test_dataclass_memory():
    validate_dataclass(new_store('dataclass-memory', clazz=Compound))


def test_dataclass():
    validate_dataclass(new_store('dataclass', clazz=Compound, index=True))


def test_invalid_store_name():
    with contextlib.suppress(NotImplementedError):
        new_store('XXXX')
