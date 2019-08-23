""" drug name, return Compound """

from bmeg import Compound, Project
from bmeg.requests import Client
from bmeg.ioutils import read_tsv

import re
import logging
import pydash

NOFINDS = {}
NOFINDS_BIOTHINGS = {}

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


def _chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]


def _decompose(name):
    """given a name, split into an array of searchable terms"""
    try:
        name = name.decode("utf-8")
    except Exception:
        pass
    name_parts = re.split('\W+', name)
    no_punct = ' '.join(name_parts).strip()
    name_parts = no_punct.split()
    pairs = [' '.join(c) for c in _chunks(name_parts, 2)]
    logging.debug('pairs {} {} >{}<'.format(pairs, len(name_parts), no_punct))
    try:
        if len(name_parts) == 1:
            logging.debug('returning [no_punct] {}'.format([no_punct]))
            return [name] + [no_punct]
        if name_parts == no_punct.split():
            logging.debug('<<< returning [name_parts] {}'.format(name_parts))
            return [name] + ['{} {}'.format(name_parts[0], name_parts[1])] + name_parts
        if [no_punct] == pairs:
            logging.debug('returning pairs + name_parts {}'.format(pairs + name_parts))
            return [name] + pairs + name_parts
    except Exception as e:
        logging.error("_decompose {}".format(name))
        logging.exception(e)
    logging.debug('returning [name] [no_punct] + pairs + name_parts {}'.format([no_punct] + pairs + name_parts))
    return [name] + [no_punct] + pairs + name_parts


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
        logging.debug('len(hits) {}'.format(len(hits)))
        if len(hits) == 0:
            logging.debug('no hit for {}'.format(name))
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


def normalize_biothings(name):
    """
     curl 'http://mychem.info/v1/query?q=chembl.molecule_synonyms.synonyms:aspirin&fields=pubchem.cid,chembl.molecule_synonyms,chembl.molecule_chembl_id,chebi.chebi_id' | jq .
    """
    compounds = []

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
    elif name.lower().startswith("cid"):
        name = name.lower().strip("cid")
        idf = 'pubchem.cid'
    elif re.match(r'^[0-9]+$', name):
        idf = 'pubchem.cid'
    else:
        idf = None

    if idf is not None:
        url = 'http://mychem.info/v1/query?q={}:{}&fields={}&size=1'.format(idf, name, fields)
        compound, score = _process_biothings_query(name, url)
        if compound is not None:
            compounds.append((compound, score))

    name_parts = _decompose(name)
    logging.debug('name parts: {}'.format(name_parts))
    for name_part in name_parts:
        logging.debug('checking {}'.format(name_part))
        if len(name_part) < NAME_PART_MIN_LEN:
            continue
        url = 'http://mychem.info/v1/query?q=chembl.pref_name:{}&fields={}&size=1'.format(name_part, fields)
        compound, score = _process_biothings_query(name_part, url)
        if compound is not None:
            compounds.append((compound, score))
        else:
            url = 'http://mychem.info/v1/query?q=chembl.molecule_synonyms.synonyms:{}&fields={}&size=1'.format(name_part, fields)
            compound, score = _process_biothings_query(name_part, url)
            if compound is not None:
                compounds.append((compound, score))
            else:
                if name_part != name:
                    NOFINDS_BIOTHINGS[name_part] = True

    if len(compounds) == 0:
        NOFINDS_BIOTHINGS[name] = True
        return None

    # sort by score and return the top hit
    compounds = sorted(compounds, key=lambda t: t[1], reverse=True)
    return compounds[0][0]


def normalize(name):
    """ given a drug name """

    if name == "N/A":
        return None
    if name in NOFINDS:
        return None

    # ensure name is a string
    try:
        name = name.decode('utf8')
    except Exception:
        pass
    name = str(name).strip()

    # do we have a better name?
    if ALIASES.get(name, None):
        logging.debug('The alias was {}'.format(ALIASES.get(name)))
    else:
        logging.debug('There was no alias for {}'.format(name))
    name = ALIASES.get(name, name)

    compound = normalize_biothings(name)
    if compound is None:
        logging.warning('normalize_drugs NOFIND {}'.format(name))
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
