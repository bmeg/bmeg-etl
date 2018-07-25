import requests
from requests.adapters import HTTPAdapter
from requests.packages.urllib3.util.retry import Retry


def Client(retries=3, backoff_factor=0.3, status_forcelist=(500, 502, 504), session=None):
    """
    Client provides a requests session that is configured for automatic retries
    for obvious failures, with exponential backoff.
    """
    # TODO include rate limiting.
    session = session or requests.Session()
    retry = Retry(
        total=retries,
        read=retries,
        connect=retries,
        backoff_factor=backoff_factor,
        status_forcelist=status_forcelist,
    )
    adapter = HTTPAdapter(max_retries=retry)
    session.mount('http://', adapter)
    session.mount('https://', adapter)
    return session
