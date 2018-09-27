from bmeg.requests import Client
import json
import logging

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
    # With a GET request, the filters parameter needs to be converted
    # from a dictionary to JSON-formatted string
    if 'filters' in params:
        params['filters'] = json.dumps(params['filters'])

    # Iterate through all the pages.
    while True:
        try:
            req = client.get(URL_BASE + endpoint, params=params)
            data = req.json()
            logging.warning(data)
            data = data['data']

            hits = data.get("hits", [])
            if len(hits) == 0:
                return

            for hit in hits:
                yield hit

            # Get the next page.
            params['from'] = data['pagination']['from'] + page_size
        except Exception as e:
            logging.warning(str(e))
            logging.warning(json.dumps(params))
            raise


def extract(data, keys):
    return {key: data.get(key) for key in keys}


def get_file(file_id, path):
    """ download a file from gdc, save in path """
    endpoint = 'data/{}'.format(file_id)
    req = client.get(URL_BASE + endpoint)
    with open(path, 'wb') as out:
        out.write(req.content)
    return path
