
"""
a set of very basic queries - simply ensure the counts of label
"""

from gripql import eq

EXPECTED_COUNTS = {
    'Individual': 34980,
    'Biosample': 74457,
    'Project': 84,
    'Aliquot': 188038,
    'DrugResponse': 197690,
    'G2PAssociation': 46573,
    'Gene': 57905,
    'MinimalAllele': 1318,
    'Allele': 4028128,
    'Phenotype': 1413,
    'Exon': 674466,
    'Protein': 107844,
    'PFAMFamily': 16712,
    'Transcript': 194318,
}


def count_label(label, V):
    """ count label template query """
    q = V.where(
        eq("_label", label)
    ).count()
    return list(
        q
    )[0]['count'], q


def test_expected_vertex_counts(V):
    """ iterate through EXPECTED_COUNTS, assert expected_count """
    errors = []
    for label in EXPECTED_COUNTS.keys():
        expected_count = EXPECTED_COUNTS[label]
        actual_count, q = count_label(label, V)
        if actual_count != expected_count:
            errors.append(
                'Expected {} {}, got {} q:{}'.format(expected_count, label, actual_count, q.query)
            )
    if len(errors) != 0:
        print(errors)
        assert False, 'Expected no errors'
