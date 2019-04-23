
"""
a set of very basic queries - simply ensure the counts of label
"""


EXPECTED_COUNTS = {
    'Case': 35596,
    'Sample': 75127,
    'Project': 196,
    'Aliquot': 195804,
    'DrugResponse': 596490,
    'G2PAssociation': 50099,
    'Gene': 63677,
    'GenomicFeature': 5941,
    'Allele': 4129937,
    'Phenotype': 3447,
    'Exon': 737019,
    'Protein': 94575,
    'PFAMFamily': 17929,
    'Transcript': 214804,
}


def count_label(label, V):
    """ count label template query """
    q = V.hasLabel(label).count()
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
