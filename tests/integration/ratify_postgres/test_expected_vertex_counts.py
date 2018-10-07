
"""
a set of very basic queries - simply ensure the counts of label
"""

EXPECTED_COUNTS = {
    'Individual': 34353,
    'Biosample': 73830,
    'Project': 65,
    'Aliquot': 187411,
    'DrugResponse': 197690,
    'G2PAssociation': 46573,
    'Gene': 26998,
    'MinimalAllele': 1318,
    'Allele': 3984706,
    'Phenotype': 1413,
    'Exon': 674466,
    'Protein': 107844,
    'PFAMFamily': 16712,
    'Transcript': 95160,
}


def count_label(label, postgres):
    """ count label template query """
    sql = "select count(*) from vertex where label = '{}'".format(label)
    return [x for x in postgres.query(sql)][0]['count'], sql


def test_expected_vertex_counts(postgres):
    """ iterate through EXPECTED_COUNTS, assert expected_count """
    errors = []
    for label in EXPECTED_COUNTS.keys():
        expected_count = EXPECTED_COUNTS[label]
        actual_count, q = count_label(label, postgres)
        if actual_count != expected_count:
            errors.append(
                'Expected {} {}, got {} q:{}'.format(expected_count, label, actual_count, q)
            )
    if len(errors) != 0:
        print(errors)
        assert False, 'Expected no errors'
