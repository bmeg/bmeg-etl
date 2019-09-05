""" drug name, return Compound """

from bmeg import Compound, Project
from bmeg.requests import Client
from bmeg.ioutils import read_tsv

import re
import logging
import pydash

NOFINDS = {}
NOFINDS_BIOTHINGS = {}
NOFINDS_PUBCHEM = {}

requests = Client('drug_enricher')

NAME_PART_MIN_LEN = 5
MIN_BIOTHINGS_SCORE = 8.5

ALIASES = {}

# get aliases
for line in read_tsv('source/drug_enricher/drug_alias.tsv'):
    if line['alias'] == 'NO-FIND':
        NOFINDS[line['name']] = True
        continue
    ALIASES[line['name']] = line['alias']


def compound_factory(name):
    """ create a stub compound for downstream normalization """
    return Compound(id=Compound.make_gid('TODO:{}'.format(name)),
                    id_source='TODO',
                    submitter_id=name,
                    project_id=Project.make_gid('Reference'))


def _process_biothings_query(name, url):
    "returns => compound: dict, score: float"

    try:
        retries = 0
        while retries < 5:
            logging.debug("_process_biothings_query: {}".format(url))
            r = requests.get(url, timeout=60)
            rsp = r.json()
            if 'hits' not in rsp:
                retries += 1
                continue
            else:
                break

        if 'hits' not in rsp:
            return None, 0.0

        hits = rsp['hits']
        logging.debug('query: {} len(hits) {}'.format(name, len(hits)))
        if len(hits) == 0:
            return None, 0.0
        # sort to get best hit
        hits = sorted(hits, key=lambda k: k['_score'], reverse=True)
        hit = hits[0]

        # hit doesn't contain the required information
        if 'pubchem' not in hit and 'chebi' not in hit and 'chembl' not in hit and 'drugbank' not in hit:
            logging.debug('no pubchem, chebi, chembl or drugbank info for {}'
                          .format(name))
            return None, 0.0

        # The higher the _score, the more relevant the document.
        if hit['_score'] < MIN_BIOTHINGS_SCORE:
            logging.debug(
                'discarded hit for {}, score too low {}'
                .format(name, hit['_score'])
            )
            return None, 0.0
        score = hit['_score']

        pubchem_id = pydash.get(hit, 'pubchem.cid', None)
        if pubchem_id:
            pubchem_id = 'CID{}'.format(pubchem_id)
        chebi_id = pydash.get(hit, 'chebi.id', None)
        chembl_id = pydash.get(hit, 'chembl.molecule_chembl_id', None)
        drugbank_id = pydash.get(hit, 'drugbank.id', None)

        ids = [pubchem_id, chebi_id, chembl_id, drugbank_id]
        if all([v is None for v in ids]):
            logging.debug('no pubchem, chebi, chembl or drugbank id for {}'
                          .format(name))
            return None, 0.0

        chembl = hit.get('chembl', {})
        synonym_fda = synonym_usan = synonym_inn = synonym_usp = None
        if 'molecule_synonyms' in chembl:
            molecule_synonyms = chembl['molecule_synonyms']
            if type(molecule_synonyms) is list:
                for molecule_synonym in molecule_synonyms:
                    if molecule_synonym['syn_type'] == 'FDA':
                        synonym_fda = molecule_synonym['synonyms'].encode('utf8')
                    if molecule_synonym['syn_type'] == 'USAN':
                        synonym_usan = molecule_synonym['synonyms'].encode('utf8')
                    if molecule_synonym['syn_type'] == 'INN':
                        synonym_inn = molecule_synonym['synonyms'].encode('utf8')
                    if molecule_synonym['syn_type'] == 'USP':
                        synonym_usp = molecule_synonym['synonyms'].encode('utf8')
            else:
                if molecule_synonyms['syn_type'] == 'FDA':
                    synonym_fda = molecule_synonyms['synonyms'].encode('utf8')
                if molecule_synonyms['syn_type'] == 'USAN':
                    synonym_usan = molecule_synonyms['synonyms'].encode('utf8')
                if molecule_synonyms['syn_type'] == 'INN':
                    synonym_inn = molecule_synonyms['synonyms'].encode('utf8')
                if molecule_synonyms['syn_type'] == 'USP':
                    synonym_usp = molecule_synonyms['synonyms'].encode('utf8')

        taxonomy = pydash.get(
            hit,
            'drugbank.taxonomy',
            None
        )
        usan_stem = pydash.get(
            hit,
            'chembl.usan_stem_definition',
            None
        )
        inchi = pydash.get(
            hit,
            'pubchem.inchi',
            pydash.get(
                hit,
                'chebi.inchi',
                pydash.get(
                    hit,
                    'chembl.inchi',
                    pydash.get(
                        hit,
                        'drugbank.inchi',
                        None
                    )
                )
            )
        )
        inchi_key = pydash.get(
            hit,
            'pubchem.inchi_key',
            pydash.get(
                hit,
                'chebi.inchi_key',
                pydash.get(
                    hit,
                    'chembl.inchi_key',
                    pydash.get(
                        hit,
                        'drugbank.inchi_key',
                        None
                    )
                )
            )
        )

        approved_countries = []
        products = pydash.get(hit, 'drugbank.products', [])
        if type(products) is list:
            for product in products:
                if product['approved']:
                    approved_countries.append(product['country'])
        else:
            product = products
            if product['approved']:
                approved_countries.append(product['country'])
        approved_countries = sorted(list(set(approved_countries)))

        id_term = pubchem_id or chembl_id or chebi_id or drugbank_id
        id_source = None
        if pubchem_id:
            id_source = "PUBCHEM"
        elif chembl_id:
            id_source = "CHEMBL"
        elif chebi_id:
            id_source = "CHEBI"
        elif drugbank_id:
            id_source = "DRUGBANK"

        compound = {
            'id': id_term,
            'id_source': id_source,
            'pubchem_id': pubchem_id,
            'chebi_id': chebi_id,
            'chembl_id': chembl_id,
            'drugbank_id': drugbank_id,
            'synonym': synonym_fda or synonym_usan or synonym_inn or synonym_usp or name,
            'inchi': inchi,
            'inchi_key': inchi_key,
            'taxonomy': taxonomy,
            'approved_countries': approved_countries,
            'usan_stem_definition': usan_stem,
            'source_url': url,
            # 'search_term': name,
        }

    except Exception as e:
        logging.exception(e)
        return None, 0.0

    return compound, score


def normalize_biothings(name, fuzzy=False):
    """
     curl 'http://mychem.info/v1/query?q=chembl.molecule_synonyms.synonyms:aspirin&fields=pubchem.cid,chembl.molecule_synonyms,chembl.molecule_chembl_id,chebi.chebi_id' | jq .
    """

    if name in NOFINDS_BIOTHINGS:
        logging.info("NOFINDS_BIOTHINGS {}".format(name))
        return None

    fields = [
        'chebi.id', 'chebi.inchi', 'chebi.inchi_key', 'chebi.name',
        'chembl.molecule_chembl_id', 'chembl.pref_name', 'chembl.inchi', 'chembl.inchi_key', 'chembl.molecule_synonyms', 'chembl.usan_stem_definition',
        'pubchem.cid', 'pubchem.inchi', 'pubchem.inchi_key',
        'drugbank.id', 'drugbank.inchi', 'drugbank.inchi_key',
        'drugbank.products.approved', 'drugbank.products.country',
        'drugbank.taxonomy.class', 'drugbank.taxonomy.direct-parent',
        'drugbank.taxonomy.kingdom', 'drugbank.taxonomy.subclass',
        'drugbank.taxonomy.superclass', 'drugbank.taxonomy.description'
    ]
    fields = ','.join(fields)

    # search for exact matches first
    # id_fields = ['chebi.id', 'chembl.molecule_chembl_id', 'pubchem.cid']
    if name.lower().startswith("chebi"):
        idf = 'chebi.id'
    elif name.lower().startswith("chembl"):
        idf = 'chembl.molecule_chembl_id'
    elif name.lower().startswith("db"):
        idf = 'drugbank.id'
    elif name.lower().startswith("cid"):
        name = re.sub('cid([\ \-_:]+)?', '', name.lower())
        idf = 'pubchem.cid'
    elif re.match(r'^[0-9]+$', name):
        idf = 'pubchem.cid'
    else:
        idf = None

    if idf is not None:
        url = 'http://mychem.info/v1/query?q={}:"{}"&fields={}&size=1'.format(idf, name, fields)
        compound, score = _process_biothings_query(name, url)
        if compound is not None:
            return compound

    # search for matches
    search_fields = [
        # names
        'chembl.pref_name',
        'chebi.name',
        'drugbank.name',
        'pharmgkb.name',
        # lists of aliases
        'chebi.brand_names',
        'chebi.synonyms'
        'pharmgkb.trade_names',
        'pharmgkb.generic_names',
        'drugbank.synonyms'
        'chembl.molecule_synonyms.synonyms',
    ]

    fuzzy_match = "~" if fuzzy else ""
    # search for a good match
    for sfield in search_fields:
        name_parts = re.split('\W+', name)
        # search for exact match in field if multiple words are present
        if len(name_parts) > 1:
            url = 'http://mychem.info/v1/query?q={}:"{}"&fields={}&size=1'.format(sfield, name, fields)
            compound, score = _process_biothings_query(name, url)
            if compound is not None:
                return compound

        # search for all name parts in field
        query_term = " AND ".join([
            "{}:{}{}".format(sfield, n, fuzzy_match) for n in name_parts
        ])
        url = 'http://mychem.info/v1/query?q={}&fields={}&size=1'.format(query_term, fields)
        compound, score = _process_biothings_query(name, url)
        if compound is not None:
            return compound

    NOFINDS_BIOTHINGS[name] = True
    return None


def search_pubchem(name):
    """
    seach pubchem and retrieve compound_id and most common synonym
    """
    if name in NOFINDS_PUBCHEM:
        logging.info("NOFINDS_PUBCHEM {}".format(name))
        return None

    search_term = name.lower()
    url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{}/cids/JSON'.format(search_term)
    if search_term.startswith("cid"):
        search_term = re.sub('cid([\ \-_:]+)?', '', search_term)
        url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{}/cids/JSON'.format(search_term)
        r = requests.get(url, timeout=60)
        rsp = r.json()
    elif search_term.startswith("sid"):
        search_term = re.sub('sid([\ \-_:]+)?', '', search_term)
        url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/sid/{}/cids/JSON'.format(search_term)

    logging.debug('search pubchem: {}'.format(url))
    r = requests.get(url, timeout=60)
    rsp = r.json()
    info = pydash.get(rsp, 'InformationList.Information.0', None)
    if info:
        rsp = {'IdentifierList': info}
    logging.debug("response: {}".format(rsp))
    cids = pydash.get(rsp, 'IdentifierList.CID', [])
    if len(cids) > 0:
        return 'CID{}'.format(cids[0])

    NOFINDS_PUBCHEM[name] = True
    return None


def normalize(name):
    """ given a drug name """

    if name == "N/A":
        return None
    if name in NOFINDS:
        logging.warning('NOFINDS {}'.format(name))
        return None

    # ensure name is a string
    try:
        name = name.decode('utf8')
    except Exception:
        pass
    name = str(name).strip()

    # do we have a better name?
    if ALIASES.get(name, None):
        logging.debug('renamed {} to {}'.format(name, ALIASES.get(name)))
    name = ALIASES.get(name, name)

    compound = normalize_biothings(name, fuzzy=False)
    if compound is None:
        cid = search_pubchem(name)
        if cid:
            compound = normalize_biothings(cid)

    if compound is None:
        compound = normalize_biothings(name, fuzzy=True)

    if compound is None:
        drug_name = re.search("([A-Za-z0-9]+)(\ +?\(.*\))?", name).group(1)
        if drug_name != name:
            compound = normalize(drug_name)

    if compound is None:
        # skip next time
        NOFINDS[name] = True
        return None

    # ensure synonym is a str, not bytes
    try:
        compound['synonym'] = compound['synonym'].decode()
    except AttributeError:
        pass

    return compound


def spell_check(name):
    """ see if google knows a better spelling, returns list of suggestions """
    url = 'http://suggestqueries.google.com/complete/search?client=firefox&q={}'.format(name)
    r = requests.get(url)
    spellings = r.json()
    suggestions = spellings[1]
    return suggestions
