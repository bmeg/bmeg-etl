from bmeg.ioutils import reader
from bmeg.util.cli import default_argument_parser
from bmeg.util.logging import default_logging
from bmeg.requests import Client
from bmeg.ioutils import read_tsv

import logging
import glob
import sys
import ujson
import gzip
from collections import defaultdict
from slugify import slugify


def open_gzip(fname):
    """Opens fname as gz."""
    return gzip.GzipFile(
        filename='',
        compresslevel=1,
        fileobj=open(fname + '.gz', 'wb'),
        mtime=0
    )


def write_page(fh, category, title, text):
    """Writes page as json."""
    page = {'category': category, 'title': title, 'text': text}
    fh.write(ujson.dumps(page).encode())
    fh.write("\n".encode())


requests = Client('drug_enricher')

REVERSE_ALIASES = defaultdict(list)
# get aliases
for line in read_tsv('source/drug_enricher/drug_alias.tsv'):
    REVERSE_ALIASES[line['alias'].lower()].append(line['name'])


TERM_ID_URL_FORMATS = {
    "INCHI": "http://rdf.ncbi.nlm.nih.gov/pubchem/inchikey/{}",
    "CHEBI": "http://purl.obolibrary.org/obo/{}",
    "CHEMBL": "https://www.ebi.ac.uk/chembl/compound_report_card/{}/",
    "PUBCHEM": "http://rdf.ncbi.nlm.nih.gov/pubchem/compound/{}",
    "DRUGBANK": "https://www.drugbank.ca/drugs/DB{}",
    "NO_ONTOLOGY": "https://www.google.com/search?q={}"
}


TERM_ID_KEY_FORMATS = {
    "INCHI": lambda k: k,
    "CHEBI": lambda k: k.replace(':', '_'),
    "CHEMBL": lambda k: k,
    "PUBCHEM": lambda k: k,
    "DRUGBANK": lambda k: k.replace('DB', ''),
    "NO_ONTOLOGY": lambda k: k,
}


def make_termId(ontology, key):
    if ontology not in TERM_ID_URL_FORMATS:
        raise Exception(f">{ontology}< not found")
    return TERM_ID_URL_FORMATS[ontology].format(TERM_ID_KEY_FORMATS[ontology](key))


def make_ontologyTitle(ontology, key):
    return f"{ontology.split(':')[-1]}-{key}"


def get_synonyms(data):
    if 'source_url' not in data:
        return []
    synonyms = set()
    hit = requests.get(data['source_url']).json()['hits'][0]
    for k in hit:
        if not type(hit[k]) == dict:
            continue
        if k == 'chebi':
            # ignore chebi.name inchi string
            continue
        for p in hit[k]:
            if 'molecule_synonyms' == p:
                if type(hit[k][p]) == dict:
                    synonyms.add(str(hit[k][p]['synonyms']))
                    continue
                for molecule_synonym in hit[k][p]:
                    synonyms.add(str(molecule_synonym['synonyms']))
                continue
            if 'name' not in p:
                continue
            synonyms.add(str(hit[k][p]))
    aliases = []
    for s in synonyms:
        if s.lower() in REVERSE_ALIASES:
            aliases.extend(REVERSE_ALIASES[s.lower()])
    if len(aliases) > 0:
        logging.debug(f"Added {aliases}")
        synonyms.update(aliases)

    return list(set([s.lower() for s in synonyms]))


def get_term(annotation):
    id = annotation['has_id']
    term = None
    for v in annotation['terms']['synonyms']:
        term = v
    if not term:
        term = annotation['terms'].get('synonym', None)
    if not term:
        term = id.split(':')[-1]
    return term


def transform(vertex_names="outputs/compound/normalized.Compound.Vertex.json*",
              output_dir="outputs",
              emitter_directory="wiki"):
    file = f'{output_dir}/{emitter_directory}/Compound.wiki.json'
    logging.info(f"writing {file}")
    fh = open_gzip(file)
    vertex_files = [filename for filename in glob.iglob(vertex_names, recursive=True)]
    ontology_keys = ['pubchem_id', 'chebi_id', 'chembl_id', 'drugbank_id', 'inchi_key']
    term_keys = ['synonym', 'submitter_id']
    for file in vertex_files:
        print(f"reading {file}")
        with reader(file) as ins:
            for line in ins:
                compound = ujson.loads(line)
                data = compound['data']
                annotation = {
                    'is_a': 'Compound',
                    'has_id': compound['gid'],
                    'belongs_to': data.get('id_source', None),
                    'ontologies': {p.replace('_id', '').replace('_key', ''): data.get(p, None) for p in ontology_keys if data.get(p, None)},
                    'terms': {p: data.get(p, None) for p in term_keys if data.get(p, None)}
                }
                annotation['terms']['synonyms'] = get_synonyms(data)
                term = get_term(annotation)
                annotation['terms']['synonyms'] = [s for s in get_synonyms(data) if not s.lower() == term.lower()]
                _synonyms = annotation['terms']['synonyms']
                synonyms = "\n".join(set([f"* synonym: [[CompoundSynonym::{slugify(s)}]]" for s in _synonyms]))

                term_ids = []
                for belongs_to, term_id in annotation['ontologies'].items():
                    belongs_to = belongs_to.upper()
                    ontology_title = make_ontologyTitle(term_id, belongs_to)
                    term_ids.append(ontology_title)
                    write_page(fh, belongs_to, ontology_title, f"""[[Category:{belongs_to.upper()}]]\n* term: [[Term::{term}]]\n* termId: [[TermId::{make_termId(belongs_to, term_id)}]]\n* compound: [[Compound::{slugify(term)}]]""")
                gid = compound['gid']

                ontology_strings = []
                for term_id in term_ids:
                    property = 'Ontology'
                    label = 'ontology'
                    k = term_id.split('-')[0]
                    if k in gid:
                        property = 'SelectedOntology'
                        label = 'selected_ontology'
                    ontology_strings.append(f"* {label}: [[{property}::{term_id}]]")
                ontologies = "\n".join(ontology_strings)
                write_page(fh, 'Compound', slugify(term), f"""[[Category:Compound]]\n* GID: [[GID::{gid}]]\n* term: [[Term::{term}]]\n{ontologies}\n{synonyms}""")
                for s in _synonyms:
                    if s.lower() == term.lower():
                        continue
                    write_page(fh, 'CompoundSynonym', slugify(s), f"""[[Category:CompoundSynonym]]\n* term: [[Term::{s}]]\n* compound: [[Compound::{slugify(term)}]]""")
    fh.close()


if __name__ == '__main__':  # pragma: no cover
    parser = default_argument_parser()
    options = parser.parse_args(sys.argv[1:])
    default_logging(options.loglevel)
    transform()
