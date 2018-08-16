""" test requests customization """

import os
import json
import contextlib
import requests
from bmeg.requests import Client

# Client will create a cache at this location
CACHE_PATH = 'outputs/test/requests_cache.sqlite'


def _remove_cache():
    """ remove cache if it exists """
    with contextlib.suppress(FileNotFoundError):
        os.remove(CACHE_PATH)


def test_Client_session():
    """ ensure cache does not get created if we supply session """
    _remove_cache()
    s = requests.Session()
    s = Client("test", session=s)
    r = s.get('http://myvariant.info/v1/variant/chr6:g.152708291G>A')
    assert json.loads(r.text), 'should deliver a response'
    assert not os.path.isfile(CACHE_PATH), 'should not create a cache'


def test_Client_cached_session():
    """ ensure cache created if we do not supply session """
    _remove_cache()
    s = Client("test")
    r = s.get('http://myvariant.info/v1/variant/chr6:g.152708291G>A')
    assert json.loads(r.text), 'should deliver a response'
    assert os.path.isfile(CACHE_PATH), 'should create a cache'

# TODO - add tests for retry
