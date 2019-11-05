from fetch_sqlite import fetch_all, sqlite_connection, features
from load_sqlite import load_table
from bmeg.requests import Client
from collections import OrderedDict


def fetch_refsnp(source_path, conn):
    """Yields refsnp annotations for each location."""
    requests = Client('pharmgkb', cache_name="ensembl_requests_cache")
    already_fetched = set()
    for annotation in fetch_all(source_path, conn, include_features=False):
        for l in features(annotation, conn):
            if l in already_fetched:
                continue
            url = f"https://grch37.rest.ensembl.org/variation/human/{l}?phenotypes=1"
            response = requests.get(url, headers={'Content-Type': 'application/json'})
            assert response.status_code in [200, 400], f"Should return 200 or 400, was {response.status_code} {url}"
            payload = response.json()
            if 'mappings' in payload:
                yield(OrderedDict({'id': l, 'data': payload}))
            else:
                print(f"WARNING: {l} returned no data {url}")
            already_fetched.add(l)


def fetch(source_path="source/pharmgkb",
          emitter_prefix=None,
          emitter_directory='pharmgkb'):
    # expects a dict
    def refsnp_cb(line):
        return line['data']
    conn = sqlite_connection(source_path)
    load_table(conn, fetch_refsnp(source_path=source_path, conn=conn), 'refsnp', refsnp_cb)


if __name__ == '__main__':  # pragma: no cover
    fetch()
