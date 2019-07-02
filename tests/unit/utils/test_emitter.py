import os
import time
import hashlib
import shutil
from bmeg.emitter import JSONEmitter
from bmeg import Compound


def test_gzip_emitter_md5():
    """Two different files, with same content, should have same hash."""
    name = 'test'
    compound = Compound(term_id='TODO:{}'.format(name), term='TODO', name=name, project_id=name)
    compound.id = Compound.make_gid(compound.term_id)
    path1 = _path(_dir(), 'test')
    emitter = JSONEmitter(path1)
    emitter.emit_vertex(compound)
    emitter.close()
    time.sleep(1)
    path2 = _path(_dir(), 'test2')
    emitter = JSONEmitter(path2)
    emitter.emit_vertex(compound)
    emitter.close()
    h1 = hashlib.md5()
    h2 = hashlib.md5()
    with open(path1 + '/Compound.Vertex.json.gz', mode='rb') as input:
        h1.update(input.read())
    with open(path2 + '/Compound.Vertex.json.gz', mode='rb') as input:
        h2.update(input.read())
    assert h1.hexdigest() == h2.hexdigest(), 'Should have the same hash'
    shutil.rmtree(path1)
    shutil.rmtree(path2)


def _dir():
    """Return the directory of this file."""
    return os.path.dirname(os.path.realpath(__file__))


def _path(*args):
    """Join args as path."""
    return os.path.join(*args)
