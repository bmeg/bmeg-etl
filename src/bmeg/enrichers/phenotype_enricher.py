""" given disease name, return Phenotype """

import logging
import re
import pydash

from bmeg import Phenotype, Project
from bmeg.requests import Client
requests = Client('phenotype_enricher')

NOFINDS = {}


def phenotype_factory(name):
    """ create a stub compound for downstream normalization """
    return Phenotype(id=Phenotype.make_gid('TODO:{}'.format(name)),
                     term_id='TODO:{}'.format(name),
                     term='TODO',
                     name=name,
                     project_id=Project.make_gid("Reference"))


def normalize_biothings(name):
    fields = ["mondo.label"]
    fields = ','.join(fields)

    # exact name in label
    url = "http://mydisease.info/v1/query?q=mondo.label:\"{}\"&fields={}&size=5".format(name, fields)
    pheno = _normalize_biothings(name, url)
    if pheno:
        return pheno

    name_parts = re.split('\W+', name)

    if len(name_parts) > 1:
        # all name parts in label
        query_term = " AND ".join(["mondo.label:{}".format(n) for n in name_parts])
        url = "http://mydisease.info/v1/query?q={}&fields={}&size=5".format(query_term, fields)
        pheno = _normalize_biothings(name, url)
        if pheno:
            return pheno

        # fuzzy match
        query_term = " AND ".join(["mondo.label:{}~".format(n) for n in name_parts])
        url = "http://mydisease.info/v1/query?q={}&fields={}&size=5".format(query_term, fields)
        pheno = _normalize_biothings(name, url)
        if pheno:
            return pheno

        # OR match
        query_term = " OR ".join(["mondo.label:{}".format(n) for n in name_parts])
        url = "http://mydisease.info/v1/query?q={}&fields={}&size=5".format(query_term, fields)
        pheno = _normalize_biothings(name, url)
        if pheno:
            return pheno

        # OR fuzzy
        # query_term = " OR ".join(["mondo.label:{}~".format(n) for n in name_parts])
        # url = "http://mydisease.info/v1/query?q={}&fields={}&size=5".format(query_term, fields)
        # pheno = _normalize_biothings(name, url)
        # if pheno:
        #     return pheno

    else:
        # fuzzy match
        url = "http://mydisease.info/v1/query?q=mondo.label:{}~&fields={}&size=5".format(name, fields)
        pheno = _normalize_biothings(name, url)
        if pheno:
            return pheno

    # exact match in any field
    url = "http://mydisease.info/v1/query?q=\"{}\"&fields={}&size=5".format(name, fields)
    pheno = _normalize_biothings(name, url)
    if pheno:
        return pheno

    # hail mary
    # url = "http://mydisease.info/v1/query?q={}&fields={}&size=5".format(name, fields)
    # pheno = _normalize_biothings(name, url)
    # if pheno:
    #     return pheno

    return None


def _normalize_biothings(name, url):
    # TODO figure out score threshold
    min_score = 7
    disease = None

    try:
        retries = 0
        while retries < 5:
            logging.debug("normalize_biothings: {}".format(url))
            r = requests.get(url, timeout=60)
            rsp = r.json()
            if 'hits' not in rsp:
                retries += 1
                continue
            else:
                break

        if 'hits' not in rsp:
            return None

        hits = rsp['hits']
        logging.debug('query: {} len(hits) {}'.format(name, len(hits)))
        if len(hits) == 0:
            return None
        # sort to get best hit
        hits = [h for h in hits if h["_id"].startswith("MONDO")]
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
        disease = {"id": hit["_id"], "term": pydash.get(hit, "mondo.label")}

    except Exception as e:
        logging.exception(e)
        return None

    print(name, disease, score)
    return disease


def normalize(name):
    name = name.strip()
    if name in NOFINDS:
        return None

    original_name = name
    name = lookup_alias(name)

    if name == "ignore":
        logging.debug('skipping {}'.format(name))
        return None

    normalized = normalize_biothings(name)
    if normalized is None:
        logging.warning('normalize_phenotype NOFIND {}'.format(original_name))
        # skip next time
        NOFINDS[original_name] = True
        return None

    return normalized


disease_alias = {}
with open('source/phenotype_enricher/disease_alias.tsv', "r") as f:
    for line in f:
        if line.startswith("#"):
            continue
        inline_list = line.rstrip().split('\t')
        disease_alias[inline_list[0].lower()] = inline_list[1]


def lookup_alias(name):
    disease = disease_alias.get(name.lower(), name)
    if disease != name:
        logging.debug('found alias for {}: {}'.format(name, disease))
    return disease.lower()
