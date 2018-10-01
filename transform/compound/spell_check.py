"""
maintains drug alias
"""

from bmeg.ioutils import reader
from bmeg.enrichers.drug_enricher import spell_check, _decompose, NAME_PART_MIN_LEN, normalize
import ujson
from nltk.metrics import edit_distance
import sys


def transform(path='outputs/compound/normalized.Compound.Vertex.json', store_path='source/drug_enricher/aliases.db'):
    suggestions = {}
    with reader(path) as ins:
        for line in ins:
            compound = ujson.loads(line)
            if 'NO_ONTOLOGY' not in compound['gid']:
                continue
            original_name = compound['data']['name']
            name_parts = _decompose(original_name)
            for name_part in name_parts:
                if len(name_part) < NAME_PART_MIN_LEN:
                    continue
                suggestion = next(iter(spell_check(name_part)), None)  # might return []
                if not suggestion:
                    continue
                distance = edit_distance(name_part.lower(), suggestion)
                if original_name.lower() == suggestion:
                    distance = 9999
                else:
                    if len(normalize(suggestion)) == 0:
                        distance = 9999
                if distance < 11:
                    suggestions[name_part] = suggestion
                else:
                    suggestions[name_part] = 'NO-FIND'

    sys.stdout.write('name{}alias{}'.format('\t', '\n'))
    for k in sorted(suggestions.keys()):
        sys.stdout.write('{}{}{}{}'.format(k.strip(), '\t', suggestions[k].strip(), '\n'))


if __name__ == '__main__':
    transform()
