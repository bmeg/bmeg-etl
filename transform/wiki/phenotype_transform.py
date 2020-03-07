from bmeg.ioutils import reader
from bmeg.util.cli import default_argument_parser
from bmeg.util.logging import default_logging
from bmeg.requests import Client
from slugify import slugify


import logging
import glob
import sys
import ujson
import gzip
from collections import defaultdict


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


requests = Client('phenotype_enricher')

# term: list of aliases
REVERSE_ALIASES = defaultdict(list)
# get aliases
with open('source/phenotype_enricher/disease_alias.tsv', "r") as f:
    for line in f:
        if line.startswith("#"):
            continue
        inline_list = line.rstrip().split('\t')
        REVERSE_ALIASES[inline_list[1].lower()].append(inline_list[0])

# 'tongue cancer': ['Tumour of tongue'],
# 'triple-receptor negative breast cancer': ['Triple Negative '
#                                            'Breast Cancer'],
# 'tuberculosis': ['TB - Tuberculosis'],


def get_mondo(term_id):
    """Returns (synonyms, ontologies) of the term."""
    synonyms = set()
    ontologies = set()
    if "mondo" in term_id.lower():
        term_id = term_id.replace('-', ':')
        url = f"https://scigraph-ontology.monarchinitiative.org/scigraph/dynamic/cliqueLeader/{term_id}.json"
        r = requests.get(url)
        assert r.status_code == 200, f"Should return a 200 {url} {r.text}"
        response = r.json()
        for n in response['nodes']:
            if not n['id'].lower() == term_id.lower():
                continue
            if 'meta' not in n:
                continue
            for s in n['meta'].get('synonym', []):
                synonyms.add(s)
            for s in n['meta'].get('http://www.geneontology.org/formats/oboInOwl#hasExactSynonym', []):
                synonyms.add(s)
            for o in n['meta'].get('http://www.geneontology.org/formats/oboInOwl#hasDbXref', []):
                ontologies.add(o)
    return (list(synonyms), list(ontologies))


def get_synonyms(annotation):
    t = annotation['term']
    synonyms = set()
    if t.lower() in REVERSE_ALIASES:
        synonyms.update(REVERSE_ALIASES[t.lower()])
    return list(synonyms)


def term_id_link(term_id):
    if "mondo" in term_id.lower():
        return f"https://monarchinitiative.org/disease/{term_id.replace('-', ':')}"
    return term_id


def transform(vertex_names="outputs/phenotype/normalized.Phenotype.Vertex.json*",
              output_dir="outputs",
              emitter_directory="wiki"):
    file = f'{output_dir}/{emitter_directory}/Phenotype.wiki.json'
    logging.info(f"writing {file}")
    fh = open_gzip(file)
    vertex_files = [filename for filename in glob.iglob(vertex_names, recursive=True)]
    for file in vertex_files:
        logging.info(f"reading {file}")
        with reader(file) as ins:
            for line in ins:
                phenotype = ujson.loads(line)
                # pprint(phenotype)
                # {'_id': 'Phenotype:MONDO:0008315',
                #  'data': {'project_id': 'Project:Reference',
                #           'submitter_id': None,
                #           'term': 'prostate cancer',
                #           'term_id': 'MONDO:0008315'},
                #  'gid': 'Phenotype:MONDO:0008315',
                #  'label': 'Phenotype'}

                data = phenotype['data']
                term_id = data.get('term_id', '')
                belongs_to = term_id.split(':')[0]

                annotation = {
                    'is_a': 'Phenotype',
                    'gid': phenotype['gid'],
                    'belongs_to': term_id.split(':')[0],
                    'ontology': term_id,
                    'term': data.get('term', None)
                }
                term = data.get('term', None)
                if not term:
                    term = data['submitter_id']
                term_id = data['term_id'].replace(':', '-')
                (mondo_synonyms, mondo_ontologies) = get_mondo(term_id)
                _synonyms = get_synonyms(annotation) + mondo_synonyms
                synonyms = "\n".join([f"* synonym: [[PhenotypeSynonym::{slugify(s)}]]" for s in _synonyms])
                gid = annotation['gid']
                ontologies = "\n".join([f"* ontology: [[Ontology::{term_id}]]" for term_id in mondo_ontologies])
                write_page(fh, 'Phenotype', slugify(term), f"""[[Category:Phenotype]]\n* GID: [[GID::{gid}]]\n* term: [[Term::{term}]]\n* selected_ontology: [[SelectedOntology::{slugify(term_id)}]]\n{ontologies}\n{synonyms}""")
                write_page(fh, belongs_to, slugify(term_id), f"""[[Category:{belongs_to}]]\n* term [[Term::{term}]]\n* termId [[TermId::{term_id_link(term_id)}]]\n* phenotype [[Phenotype::{slugify(term)}]]""")
                for s in _synonyms:
                    write_page(fh, 'PhenotypeSynonym', slugify(s), f"""[[Category:PhenotypeSynonym]]\n* term: [[Term::{s}]]\n* phenotype: [[Phenotype::{slugify(term)}]]""")
    fh.close()


if __name__ == '__main__':  # pragma: no cover
    parser = default_argument_parser()
    options = parser.parse_args(sys.argv[1:])
    default_logging(options.loglevel)
    transform()
