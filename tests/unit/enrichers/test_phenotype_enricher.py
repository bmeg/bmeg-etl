import bmeg.enrichers.phenotype_enricher as phenotype_enricher
import logging
from pprint import pprint

# maintain fresh cache for all calls
from uncacher import uncache
uncache(phenotype_enricher.requests)


NAMES = """
BOCA
cancer
Hodgkin_lymphoma
transitional_cell_carcinoma
""".strip().split("\n")

EXPECTED = [
    {'ontology_term': 'MONDO:0002129', 'label': 'bone cancer', 'source': 'https://monarchinitiative.org/disease/MONDO:0002129', 'provenance': 'https://api.monarchinitiative.org/api/search/entity/bone cancer?category=disease&prefix=MONDO&start=0&rows=5'},
    {'ontology_term': 'MONDO:0004992', 'label': 'cancer', 'source': 'https://monarchinitiative.org/disease/MONDO:0004992', 'provenance': 'https://api.monarchinitiative.org/api/search/entity/cancer?category=disease&prefix=MONDO&start=0&rows=5'},
    {'ontology_term': 'MONDO:0004952', 'label': 'Hodgkins lymphoma', 'source': 'https://monarchinitiative.org/disease/MONDO:0004952', 'provenance': 'https://api.monarchinitiative.org/api/search/entity/hodgkin_lymphoma?category=disease&prefix=MONDO&start=0&rows=5'},
    # http://www.ontobee.org/ontology/MONDO?iri=http://purl.obolibrary.org/obo/MONDO_0002828
    {'ontology_term': 'MONDO:0002828', 'label': 'Bartholin gland transitional cell carcinoma', 'source': 'https://monarchinitiative.org/disease/MONDO:0002828', 'provenance': 'https://api.monarchinitiative.org/api/search/entity/transitional_cell_carcinoma?category=disease&prefix=MONDO&start=0&rows=5'}
]


def test_simple(caplog):
    """ straightforward """
    caplog.set_level(logging.DEBUG)
    for name, expected in zip(NAMES, EXPECTED):
        actual = phenotype_enricher.normalize(name)
        try:
            assert actual, 'Should return value for {}'.format(name)
            assert len(actual) == len(expected), 'Should return same number of hits for {}'.format(name)
            for f in expected:
                assert f in actual, '{} not found in actual'.format(f)
                assert actual[f] == expected[f], '{} mismatch'.format(f)
        except Exception as e:
            print(str(e))
            pprint(actual)
            pprint(expected)
            assert False, 'Failed'
