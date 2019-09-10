""" given disease name, return Phenotype """

import logging
import os
import pydash
import unicodedata
import urllib
import re


from bmeg import Phenotype, Project
from bmeg.requests import Client
requests = Client('phenotype_enricher')

NOFINDS = {}
BIOONTOLOGY_NOFINDS = {}

API_KEY = os.environ.get('BIOONTOLOGY_API_KEY')


disease_alias = {}
with open('source/phenotype_enricher/disease_alias.tsv', "r") as f:
    for line in f:
        if line.startswith("#"):
            continue
        inline_list = line.rstrip().split('\t')
        disease_alias[inline_list[0].lower()] = inline_list[1]


def phenotype_factory(name):
    """ create a stub compound for downstream normalization """
    return Phenotype(id=Phenotype.make_gid('TODO:{}'.format(name)),
                     term_id='TODO:{}'.format(name),
                     term='TODO',
                     name=name,
                     project_id=Project.make_gid("Reference"))


def normalize_bioontology(name):
    """ call bioontology & retrieve """
    if not API_KEY:
        raise ValueError('Please set BIOONTOLOGY_API_KEY in environment')

    if name in BIOONTOLOGY_NOFINDS:
        logging.info('{} in disease_normalizer.BIOONTOLOGY_NOFINDS'
                     .format(name))
        return []
    quoted_name = urllib.parse.quote_plus(name)
    url = 'http://data.bioontology.org/search?q={}&apikey={}'.format(quoted_name, API_KEY)  # NOQA
    r = requests.get(url, timeout=20)
    r.raise_for_status()
    response = r.json()
    terms = []
    if 'collection' in response and len(response['collection']) > 0:
        collection = response['collection'][0]
        parts = collection['@id'].split('/')
        ontology = parts[-2]
        id = parts[-1]
        if ontology == 'obo':
            (ontology, id) = id.split('_')
        term = {'ontology_term': '{}:{}'.format(ontology, id),
                'label': name,
                'source': ontology,
                'provenance': url}
        terms.append(term)
        family = get_family(term['ontology_term'])
        if family:
            term['family'] = family
    else:
        BIOONTOLOGY_NOFINDS[name] = True
    return terms


def normalize_ebi(name):
    """ call ebi & retrieve """
    if name in NOFINDS:
        # logging.info('{} in disease_normalizer.NOFINDS'.format(name))
        return []
    name = urllib.parse.quote_plus(name)
    url = 'https://www.ebi.ac.uk/ols/api/search?q={}&groupField=iri&exact=on&start=0&ontology=mondo'.format(name)  # NOQA
    # .response
    """
    {
      "numFound": 1,
      "start": 0,
      "docs": [
        {
          "id": "doid:http://purl.obolibrary.org/obo/DOID_1909",
          "iri": "http://purl.obolibrary.org/obo/DOID_1909",
          "short_form": "DOID_1909",
          "obo_id": "DOID:1909",
          "label": "melanoma",
          "description": [
            "A cell type cancer that has_material_basis_in abnormally proliferating cells derives_from melanocytes which are found in skin, the bowel and the eye."
          ],
          "ontology_name": "doid",
          "ontology_prefix": "DOID",
          "type": "class",
          "is_defining_ontology": true
        }
      ]
    }
    """  # NOQA

    r = requests.get(url, timeout=20)
    r.raise_for_status()
    rsp = r.json()
    if 'response' not in rsp:
        logging.info('{} in disease_normalizer.NOFINDS'.format(name))
        NOFINDS.append(name)
        return []
    response = rsp['response']
    numFound = response['numFound']
    if numFound == 0:
        logging.info('{} in disease_normalizer.NOFINDS'.format(name))
        NOFINDS.append(name)
        return []
    doc = response['docs'][0]
    # # check whether returned info is actually DOID or some other response
    # # since we only want DOID entries
    # if doc['obo_id'][:2] != 'DO':
    #     logging.info('{} in disease_normalizer.NOFINDS'.format(name))
    #     NOFINDS.append(name)
    #     return []
    term = {'ontology_term': doc['obo_id'].encode('utf8'),
            'label': doc['label'].encode('utf8'),
            'source': doc['iri'].encode('utf8'),
            'provenance': url
            }
    family = get_family(doc['obo_id'])
    if family:
        term['family'] = family

    return [term]


def get_family(ontology_id):
    """
    get the hierarchy
    """
    if not API_KEY:
        raise ValueError('Please set BIOONTOLOGY_API_KEY in environment')

    url = r = None
    try:
        if ontology_id.startswith('DOID'):
            url = 'http://disease-ontology.org/query_tree?search=True&node={}'.format(ontology_id)  # NOQA
            r = requests.get(url, timeout=20)
            if not r.status_code == 500:
                rsp = r.json()
                return get_hierarchy_family(get_hierarchy(rsp[0], []))['text']

        (ontology, k) = ontology_id.split(':')
        if ontology not in ['SNOMEDCT', 'DOID', 'RCD', 'MESH',
                            'OMIM', 'MEDDRA', 'MEDLINEPLUS', 'RCD']:
            k = ontology_id
        url = 'http://data.bioontology.org/ontologies/{}/classes/{}/ancestors?apikey={}'.format(ontology, k, API_KEY)  # NOQA
        r = requests.get(url, timeout=20)
        if r.status_code == 200:
            classes = r.json()
            if len(classes) > 0:
                for clazz in classes:
                    if '_' not in clazz['prefLabel']:
                        if re.match('^C.[0-9]*', clazz['prefLabel']):
                            continue
                        return clazz['prefLabel']
        return None
    except Exception as e:
        logging.exception(e)
        logging.error('get_family {} {} {}'.format(url, r, e))
        return None


def get_hierarchy(node, hierarchy_list):
    if (
        node.get('iconCls', '') == 'search-select-icon' or node.get('expanded', False)
    ):
        hierarchy_list.append({'text': node['text'], 'id': node['id']})
        for n in node.get('children', []):
            get_hierarchy(n, hierarchy_list)
    return hierarchy_list


def print_hierarchy(hierarchy_list, indent=0):
    for node in hierarchy_list:
        print(''.ljust(indent), node['text'], node['id'])
        indent = indent + 2


def get_hierarchy_family(_a):
    midpoint = int(len(_a) / 2) + 1
    return _a[midpoint]


def normalize_monarch(name):
    try:
        name = unicodedata.normalize('NFD', name)\
                          .encode('ascii', 'ignore')\
                          .decode()\
                          .replace("'", "")
    except Exception:
        logging.warning("failed to normalize string")
    min_score = 35
    size = 5
    url_parts = [
        'https://api.monarchinitiative.org/api/search/entity/{}?'.format(name.replace('/', ' ')),
        'category=disease&',
        'prefix=MONDO&',
        'start=0&',
        'rows={}'.format(size)
    ]
    url = ''.join(url_parts)
    logging.debug('normalize_monarch: {}'.format(url))
    try:
        r = requests.get(url, timeout=60)
        rsp = r.json()
        if len(rsp.get('docs', [])) < 1:
            return None
        # sort to get best hit
        hits = rsp['docs']
        # hits = sorted(hits, key=lambda k: k['score'], reverse=True)
        hit = hits[0]
        # check the scores
        if hits[0]['score'] < min_score:
            logging.debug(
                'discarded hit for {}, score too low {}'
                .format(name, hits[0]['score'])
            )
            return None
        # apply some basic logic if there are several equally scored hits
        if len(hits) > 1 and hits[0]['score'] == hits[1]['score']:
            logging.debug(
                'multiple hits with same score for {}'
                .format(name)
            )
            found = False
            for h in [x for x in hits if x['score'] == hits[0]['score']]:
                if name == h['label'][0]:
                    hit = h
                    found = True
            if not found:
                for h in [x for x in hits if x['score'] == hits[0]['score']]:
                    if name in h['label'][0]:
                        hit = h
                        found = True
            if not found:
                logging.warning(
                    'ambiguous hits with same score for {}; using first hit'
                    .format(name)
                )

        disease = {'ontology_term': hit['id'],
                   'label': pydash.get(hit, 'label.0', name),
                   'source': 'https://monarchinitiative.org/disease/{}'.format(hit['id']),
                   'provenance': url}
        # print(name, disease.get('ontology_term', 'NONE'), disease.get('label', 'NONE'), hit['score'], sep='\t')
        return disease

    except Exception as e:
        logging.exception(e)

    return None


def normalize_biothings(name):
    fields = ['mondo.label']
    fields = ','.join(fields)

    # exact name in label
    url = 'http://mydisease.info/v1/query?q=mondo.label:"{}"&fields={}&size=5'.format(name, fields)
    pheno = _normalize_biothings(name, url)
    if pheno:
        return pheno

    name_parts = re.split('\W+', name)

    if len(name_parts) > 1:
        # all name parts in label
        query_term = ' AND '.join(['mondo.label:{}'.format(n) for n in name_parts])
        url = 'http://mydisease.info/v1/query?q={}&fields={}&size=5'.format(query_term, fields)
        pheno = _normalize_biothings(name, url)
        if pheno:
            return pheno

        # fuzzy match
        query_term = ' AND '.join(['mondo.label:{}~'.format(n) for n in name_parts])
        url = 'http://mydisease.info/v1/query?q={}&fields={}&size=5'.format(query_term, fields)
        pheno = _normalize_biothings(name, url)
        if pheno:
            return pheno

        # OR match
        query_term = ' OR '.join(['mondo.label:{}'.format(n) for n in name_parts])
        url = 'http://mydisease.info/v1/query?q={}&fields={}&size=5'.format(query_term, fields)
        pheno = _normalize_biothings(name, url)
        if pheno:
            return pheno

        # OR fuzzy
        # query_term = ' OR '.join(['mondo.label:{}~'.format(n) for n in name_parts])
        # url = 'http://mydisease.info/v1/query?q={}&fields={}&size=5'.format(query_term, fields)
        # pheno = _normalize_biothings(name, url)
        # if pheno:
        #     return pheno

    else:
        # fuzzy match
        url = 'http://mydisease.info/v1/query?q=mondo.label:{}~&fields={}&size=5'.format(name, fields)
        pheno = _normalize_biothings(name, url)
        if pheno:
            return pheno

    # exact match in any field
    url = 'http://mydisease.info/v1/query?q="{}"&fields={}&size=5'.format(name, fields)
    pheno = _normalize_biothings(name, url)
    if pheno:
        return pheno

    # hail mary
    # url = 'http://mydisease.info/v1/query?q={}&fields={}&size=5'.format(name, fields)
    # pheno = _normalize_biothings(name, url)
    # if pheno:
    #     return pheno

    return None


def _normalize_biothings(name, url):
    # TODO figure out score threshold
    min_score = 7.0
    disease = None

    try:
        logging.debug('normalize_biothings: {}'.format(url))
        r = requests.get(url, timeout=60)
        rsp = r.json()
        if 'hits' not in rsp:
            return None

        hits = rsp['hits']
        logging.debug('query: {} len(hits) {}'.format(name, len(hits)))
        if len(hits) == 0:
            return None
        # sort to get best hit
        hits = [h for h in hits if h['_id'].startswith('MONDO')]
        hits = sorted(hits, key=lambda k: k['_score'], reverse=True)
        hit = hits[0]
        # The higher the _score, the more relevant the document.
        if hit['_score'] < min_score:
            logging.debug(
                'discarded hit for {}, score too low {}'
                .format(name, hit['_score'])
            )
            return None
        score = hit['_score']
        disease = {'ontology_term': hit['_id'],
                   'label': pydash.get(hit, 'mondo.label'),
                   'source': None,
                   'provenance': url}

    except Exception as e:
        logging.exception(e)
        return None

    print(name, disease.get('label', 'NONE'), score, sep='\t')
    return disease


def normalize(name):
    name = name.strip()
    if name in NOFINDS:
        return None

    original_name = name
    name = lookup_alias(name)

    if name == 'ignore':
        logging.debug('skipping {}'.format(name))
        return None

    normalized = normalize_monarch(name)
    if normalized is None:
        logging.warning('normalize_phenotype NOFIND {}'.format(original_name))
        # skip next time
        NOFINDS[original_name] = True
        return None

    return normalized


def lookup_alias(name):
    disease = disease_alias.get(name.lower(), name)
    if disease != name:
        logging.debug('found alias for {}: {}'.format(name, disease))
    return disease.lower()
