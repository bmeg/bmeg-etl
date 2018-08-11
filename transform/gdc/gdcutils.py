from bmeg.requests import Client

URL_BASE = "https://api.gdc.cancer.gov/"
client = Client("gdc")


def query_gdc(endpoint, params):
    """
    query_gdc makes a query to the GDC API while handling common issues
    like pagination, retries, etc.

    The return value is an iterator.
    """
    # Copy input params to avoid modification.
    params = dict(params)
    page_size = 100
    params['size'] = page_size

    # Iterate through all the pages.
    while True:
        req = client.get(URL_BASE + endpoint, params=params)
        data = req.json()['data']

        hits = data.get("hits", [])
        if len(hits) == 0:
            return

        for hit in hits:
            yield hit

        # Get the next page.
        params['from'] = data['pagination']['from'] + page_size


def extract(data, keys):
    return {key: data.get(key) for key in keys}
