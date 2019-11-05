from load_sqlite import load_table
from bmeg.requests import Client
from fetch_sqlite import sqlite_connection


def fetch(source_path="source/pharmgkb",
          emitter_prefix=None,
          emitter_directory='pharmgkb'):

    #
    requests = Client('pharmgkb', cache_name="pharmgkb_variant_requests_cache", status_forcelist=(500, 502, 504, 429), method_whitelist=['GET', 'POST'])

    conn = sqlite_connection(source_path)

    def fetch_variant(id):
        url = f"https://api.pharmgkb.org/v1/data/variant/{id}?view=max"
        response = requests.get(url)
        assert response.status_code == 200, f"Should have 200 status {response.status_code} {url}"
        data = response.json()
        assert 'data' in data, 'Response missing "data"'
        assert 'alleles' in data['data'], 'Response missing "alleles"'
        return [a['location']['displayName'] for a in data['data']['alleles']]

    def fetch_variants():
        """Yields refsnp annotations for each location."""
        start = 0
        total = None
        while True:
            query = {"query": None, "objCls": "variant", "hasData": None, "sort": True, "size": 100, "from": start, "prefix": ""}
            url = "https://api.pharmgkb.org/v1/site/browse/variant"
            print(start, 'of', total)
            response = requests.post(url, json=query)
            assert response.status_code == 200, f"Should have 200 status {response.status_code} {url}"
            data = response.json()
            assert 'data' in data, 'Response missing "data"'
            assert 'hits' in data['data'], 'Response missing "hits"'
            if not total:
                total = data['data']['total']
            c = 0
            if len(data['data']['hits']) == 0:
                break
            for hit in data['data']['hits']:
                yield {'name': hit['name'], 'id': hit['id'], 'variant': hit}
                c += 1
            start += c

    load_table(conn, fetch_variants(), 'variants')


if __name__ == '__main__':  # pragma: no cover
    fetch()
