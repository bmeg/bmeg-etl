import os

from requests.adapters import HTTPAdapter
from requests.packages.urllib3.util.retry import Retry
import requests_cache

from bmeg.utils import ensure_directory


def Client(prefix, cache_name="requests_cache", retries=5, backoff_factor=0.3,
           status_forcelist=(500, 502, 504), session=None, method_whitelist=Retry.DEFAULT_METHOD_WHITELIST):
    """
    Client provides a requests session that is configured with automatic
    retries, caching, and more.
    """
    ensure_directory("outputs", prefix)
    cache_path = os.path.join("outputs", prefix, cache_name)
    # TODO include rate limiting.
    session = session or requests_cache.CachedSession(cache_path)
    retry = Retry(
        total=retries,
        read=retries,
        connect=retries,
        backoff_factor=backoff_factor,
        status_forcelist=status_forcelist,
        method_whitelist=method_whitelist
    )
    adapter = HTTPAdapter(max_retries=retry)
    session.mount('http://', adapter)
    session.mount('https://', adapter)
    return session
