from bmeg.ioutils import reader
from bmeg.enrichers.drug_enricher import get_pubchem_autocomplete_suggestion
import ujson
import sys


def transform(path='outputs/compound/normalized.Compound.Vertex.json.gz'):
    """
    used to create source/drug_enricher/drug_alias.tsv
    """
    sys.stdout.write('name{}alias{}'.format('\t', '\n'))
    with reader(path) as ins:
        for line in ins:
            compound = ujson.loads(line)
            if 'NO_ONTOLOGY' not in compound['gid']:
                continue
            original_name = compound['data']['submitter_id'].strip()
            suggestion = get_pubchem_autocomplete_suggestion(original_name)
            if not suggestion:
                continue

            sys.stdout.write('{}{}{}{}'.format(original_name.strip(), '\t', suggestion.strip(), '\n'))


if __name__ == '__main__':
    transform()
