
import bmeg.enrichers.drug_enricher as drug_enricher
import logging
from pprint import pprint

# maintain fresh cache for all calls
from uncacher import uncache
uncache(drug_enricher.requests)


NAMES = """
asprin
5-Fluorouracil
(5Z)-7-Oxozeaenol
A-443654
A-770041
Afatinib (1)
Afatinib (2)
AICA Ribonucleotide
AKT inhibitor VIII (1)
Ecotrin
Tomaxifen
Tamoxiten
Abagovomag
Zometa
""".strip().split("\n")

EXPECTED = [
    [{'ontology_term': 'CID2244', 'source': 'http://rdf.ncbi.nlm.nih.gov/pubchem/compound', 'synonym': 'aspirin'}],
    [{'ontology_term': 'CID3385', 'synonym': 'FLUOROURACIL', 'source': 'http://rdf.ncbi.nlm.nih.gov/pubchem/compound', 'toxicity': 'LD<sub>50</sub>=230mg/kg (orally in mice)', 'taxonomy': {'class': 'Diazines', 'direct-parent': 'Halopyrimidines', 'kingdom': 'Organic compounds', 'subclass': 'Pyrimidines and pyrimidine derivatives', 'superclass': 'Organoheterocyclic compounds'}, 'approved_countries': ['Canada', 'US'], 'usan_stem': 'uracil derivatives used as thyroid antagonists and as antineoplastics'}],
    [{'ontology_term': 'CID45485157', 'source': 'http://rdf.ncbi.nlm.nih.gov/pubchem/compound', 'synonym': '(5Z)-7-Oxozeaenol'}, {'ontology_term': 'CID9863776', 'source': 'http://rdf.ncbi.nlm.nih.gov/pubchem/compound', 'synonym': 'Oxozeaenol'}],
    [{'ontology_term': 'CID10172943', 'synonym': 'A-443654', 'source': 'http://rdf.ncbi.nlm.nih.gov/pubchem/compound', 'taxonomy': {'class': 'Indoles and derivatives', 'direct-parent': '3-alkylindoles', 'kingdom': 'Organic compounds', 'subclass': 'Indoles', 'superclass': 'Organoheterocyclic compounds'}}],
    [{'ontology_term': 'CID9549184', 'source': 'http://rdf.ncbi.nlm.nih.gov/pubchem/compound', 'synonym': 'A-770041'}],
    [{'ontology_term': 'CID10184653', 'synonym': 'AFATINIB', 'source': 'http://rdf.ncbi.nlm.nih.gov/pubchem/compound', 'toxicity': 'Most common adverse reactions (≥20%) are diarrhea, rash/dermatitis, acneiform, stomatitis, paronychia, dry skin, decreased appetite, pruritus [FDA Label].\r\n\r\nConversely, overdose in 2 healthy adolescents involving the ingestion of 360 mg each of afatinib (as part of a mixed drug ingestion) was associated with adverse events of nausea, vomiting, asthenia, dizziness, headache, abdominal pain and elevated amylase (< 1.5 times ULN) [L2937]. Both individuals recovered from these adverse events [L2937].', 'taxonomy': {'class': 'Diazanaphthalenes', 'direct-parent': 'Quinazolinamines', 'kingdom': 'Organic compounds', 'subclass': 'Benzodiazines', 'superclass': 'Organoheterocyclic compounds'}, 'approved_countries': ['Canada', 'US'], 'usan_stem': 'tyrosine kinase inhibitors'}, ],
    [{'ontology_term': 'CID10184653', 'synonym': 'AFATINIB', 'source': 'http://rdf.ncbi.nlm.nih.gov/pubchem/compound', 'toxicity': 'Most common adverse reactions (≥20%) are diarrhea, rash/dermatitis, acneiform, stomatitis, paronychia, dry skin, decreased appetite, pruritus [FDA Label].\r\n\r\nConversely, overdose in 2 healthy adolescents involving the ingestion of 360 mg each of afatinib (as part of a mixed drug ingestion) was associated with adverse events of nausea, vomiting, asthenia, dizziness, headache, abdominal pain and elevated amylase (< 1.5 times ULN) [L2937]. Both individuals recovered from these adverse events [L2937].', 'taxonomy': {'class': 'Diazanaphthalenes', 'direct-parent': 'Quinazolinamines', 'kingdom': 'Organic compounds', 'subclass': 'Benzodiazines', 'superclass': 'Organoheterocyclic compounds'}, 'approved_countries': ['Canada', 'US'], 'usan_stem': 'tyrosine kinase inhibitors'}, ],
    [{'ontology_term': 'CID65110', 'source': 'http://rdf.ncbi.nlm.nih.gov/pubchem/compound', 'synonym': 'AICA ribonucleotide'}],
    [{'ontology_term': 'SID319552803', 'source': 'http://rdf.ncbi.nlm.nih.gov/pubchem/substance', 'synonym': 'akt inhibitor'}],
    [{'ontology_term': 'CID2244', 'source': 'http://rdf.ncbi.nlm.nih.gov/pubchem/compound', 'synonym': 'aspirin'}],
    [{'approved_countries': ['Canada', 'US'], 'ontology_term': 'CID2733526', 'source': 'http://rdf.ncbi.nlm.nih.gov/pubchem/compound', 'synonym': 'Tamoxifen', 'taxonomy': {'class': 'Stilbenes', 'direct-parent': 'Stilbenes', 'kingdom': 'Organic compounds', 'superclass': 'Phenylpropanoids and polyketides'}, 'toxicity': 'Signs observed at the highest doses following studies to ' 'determine LD<sub>50</sub> in animals were respiratory ' 'difficulties and convulsions.'}],
    [{'approved_countries': ['Canada', 'US'], 'ontology_term': 'CID2733526', 'source': 'http://rdf.ncbi.nlm.nih.gov/pubchem/compound', 'synonym': 'Tamoxifen', 'taxonomy': {'class': 'Stilbenes', 'direct-parent': 'Stilbenes', 'kingdom': 'Organic compounds', 'superclass': 'Phenylpropanoids and polyketides'}, 'toxicity': 'Signs observed at the highest doses following studies to ' 'determine LD<sub>50</sub> in animals were respiratory ' 'difficulties and convulsions.'}],
    [{'ontology_term': 'CHEMBL1742981', 'source': 'http://rdf.ebi.ac.uk/terms/chembl', 'synonym': 'Abagovomab', 'usan_stem': 'monoclonal antibodies'}],
    [{'approved_countries': ['Canada', 'EU', 'US'], 'ontology_term': 'CID68740', 'source': 'http://rdf.ncbi.nlm.nih.gov/pubchem/compound', 'synonym': 'ZOLEDRONIC ACID', 'taxonomy': {'class': 'Organic phosphonic acids and derivatives', 'direct-parent': 'Bisphosphonates', 'kingdom': 'Organic compounds', 'subclass': 'Bisphosphonates', 'superclass': 'Organic acids and derivatives'}, 'toxicity': 'There is no experience of acute overdose. Two patients received ' 'zoledronate (32 mg) over 5 minutes in clinical trials. Neither ' 'patient experienced any clinical or laboratory toxicity. ' 'Overdosage may cause clinically significant hypocalcemia, ' 'hypophosphatemia, and hypomagnesemia.', 'usan_stem': 'calcium metabolism regulators'}],
]


def test_simple(caplog):
    """ straightforward """
    caplog.set_level(logging.DEBUG)

    for name, expected in zip(NAMES, EXPECTED):
        actual = drug_enricher.normalize(name)
        try:
            assert actual, 'Should return value for {}'.format(name)
            assert len(actual) == len(expected), 'Should return same number of hits for {}'.format(name)
            assert actual == expected
        except Exception as e:
            print(str(e))
            pprint(actual)
            pprint(expected)
            assert False, 'Failed'


def test_alias():
    assert 'Tomaxifen' in drug_enricher.ALIASES, 'we should have an alias for Tomaxifen'
    assert drug_enricher.ALIASES['Tomaxifen'].lower() == 'tamoxifen', 'the alias for Tomaxifen should be tamoxifen'
    assert 'Tamoxiten' in drug_enricher.ALIASES, 'we should have an alias for Tamoxiten'
    assert drug_enricher.ALIASES['Tamoxiten'].lower() == 'tamoxifen', 'the alias for Tamoxiten should be tamoxifen'


def test_spell_check():
    """ """
    assert drug_enricher.spell_check('Tamoxiten')[0] == 'tamoxifen'
