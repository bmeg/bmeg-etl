"""Deletes url from cache before GET."""

import os
import sys

_session = None
_original_get = None


# create a new get method
def _uncache_before_get(*k, **kw):
    """If cache is enabled, delete url from cache before call."""
    global _original_get, _session
    url = kw.get('url', None)
    if not url and len(k) > 0:
        url = k[0]
    if url:
        cache = None
        try:
            cache = _session.cache
        except Exception as e:
            print(e)
        if cache:
            try:
                cache.delete_url(url)
            except Exception as e:
                print('>uncacher:could not remove from cache {}'.format(str(e)), file=sys.stderr)
    return _original_get(*k, **kw)


def uncache(session):
    """ """
    global _original_get, _session
    if 'USE_CACHE' not in os.environ:
        # save the original get method
        # replace original method with our monkey patch
        _session = session
        _original_get = session.get
        session.get = _uncache_before_get
        print('>uncacher:disabled cache, set USE_CACHE environmental variable to enable.', file=sys.stderr)
