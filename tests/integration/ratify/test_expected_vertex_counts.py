
"""
a set of very basic queries - simply ensure the counts of label
"""

from gripql import eq

EXPECTED_COUNTS = {
    'Individual': 32555,
    'Biosample': 58232,
    'Project': 40,
    'Aliquot': 170767,
    'DrugResponse': 197690,
    'G2PAssociation': 43369,
    'Gene': 26998,
    'MinimalAllele': 1142,
    'Allele': 3987487,
    'Phenotype': 1357,
    'Exon': 674466,
    'Protein': 107844,
    'PFAMFamily': 16712,
    'Transcript': 95160,
}


def count_label(label, V):
    """ count label template query """
    return list(
        V.where(
            eq("_label", label)
        ).count()
    )[0]['count']


def test_expected_counts(V):
    """ iterate through EXPECTED_COUNTS, assert expected_count """
    errors = []
    for label in EXPECTED_COUNTS.keys():
        expected_count = EXPECTED_COUNTS[label]
        actual_count = count_label(label, V)
        if actual_count != expected_count:
            errors.append(
                'Expected {} {}, got {}'.format(expected_count, label, actual_count)
            )
    assert len(errors) == 0, errors
