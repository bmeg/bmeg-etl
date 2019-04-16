""" given disease name, return Phenotype """

from bmeg import Phenotype
import urllib.parse
import logging
import re
import os

from bmeg.requests import Client
requests = Client('phenotype_enricher')


NOFINDS = []
BIOONTOLOGY_NOFINDS = []

API_KEY = os.environ.get('BIOONTOLOGY_API_KEY')
if not API_KEY:
    raise ValueError('Please set BIOONTOLOGY_API_KEY in environment')


def phenotype_factory(name):
    """ create a stub compound for downstream normalization """
    return Phenotype(term_id='TODO:{}'.format(name), term='TODO', name=name)


disease_alias = {}
with open('source/phenotype_enricher/disease_alias.tsv', "r") as f:
    for line in f:
        if line.startswith("#"):
            continue
        inline_list = line.rstrip().split('\t')
        disease_alias[inline_list[0].lower()] = inline_list[1]


def normalize_bioontology(name):
    """ call bioontology & retrieve """
    if name in BIOONTOLOGY_NOFINDS:
        logging.info('{} in disease_normalizer.BIOONTOLOGY_NOFINDS'
                     .format(name))
        return []
    quoted_name = urllib.parse.quote_plus(name)
    url = 'http://data.bioontology.org/search?q={}&apikey={}'.format(quoted_name, API_KEY)  # NOQA
    r = requests.get(url, timeout=20)
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
        BIOONTOLOGY_NOFINDS.append(name)
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
    # get the hierarchy
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


def normalize(name):
    try:
        diseases = []
        original_name = name
        name = project_lookup(name)
        if name:
            # find in ebi
            normalized_diseases = normalize_ebi(name)
            if len(normalized_diseases) > 0:
                diseases = diseases + normalized_diseases
            else:
                names = re.split("[\,;_]+", name)
                names.insert(0, ' '.join(names))
                for name_part in names:
                    name_part = project_lookup(name_part)
                    logging.debug("name_part {}".format(name_part))
                    normalized_diseases = normalize_ebi(name_part)
                    diseases = diseases + normalized_diseases
            if len(diseases) == 0:
                diseases = normalize_bioontology(name)
        # add the original name back
        for d in diseases:
            d['name'] = original_name
        # dedupe
        ontology_terms = {}
        for d in diseases:
            for f in d.keys():
                try:
                    d[f] = d[f].decode()
                except AttributeError:
                    pass
            if d['ontology_term'] not in ontology_terms:
                ontology_terms[d['ontology_term']] = d
        return diseases
    except Exception as e:
        logging.exception(e)
        logging.warning("Could not normalize {}".format(name))
        return []


def project_lookup(name):
    disease = disease_alias.get(name.lower())
    if not disease == name and disease:
        logging.debug('renamed {} to {}'.format(name, disease))
    if not disease:
        disease = name
    return disease
