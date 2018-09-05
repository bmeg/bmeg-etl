
"""
a set of very basic queries - simply ensure the counts of label->label
"""


from gripql import eq

EXPECTED_COUNTS = [
    {'_from': 'Biosample', 'to': 'Individual', 'expected_count': 57186},
    {'_from': 'Individual', 'to': 'Project', 'expected_count': 32555},
    {'_from': 'Aliquot', 'to': 'Biosample', 'expected_count': 170767},
    {'_from': 'DrugResponse', 'to': 'Biosample', 'expected_count': 169865},
    {'_from': 'Protein', 'to': 'PFAMFamily', 'expected_count': 108729},
    {'_from': 'Protein', 'to': 'Transcript', 'expected_count': 73439},
]


def count_traversal(_from, to, expected_count, V):
    """ count traversal template query """
    actual_count = list(
        V.where(
            eq("_label", _from)
        ).out().where(eq("_label", to)).count()
    )[0]['count']
    if actual_count != expected_count:
        return 'Expected from:{} to:{} expected:{} actual: {}'.format(
            _from, to, expected_count, actual_count
        )


def test_expected_counts(V):
    """ iterate through EXPECTED_COUNTS, assert expected_count """
    errors = []
    for traversal in EXPECTED_COUNTS:
        error_msg = count_traversal(**traversal, V=V)
        if error_msg:
            errors.append(error_msg)
    assert len(errors) == 0, errors


def test_expected_exon_transcript(V):
    """ subset only for one chromosome """
    q = (
        V
        .where(eq("_label", 'Gene'))
        .where(eq("chromosome", '22'))
        .in_().where(eq("_label", 'Transcript'))
        .count()
    )
    actual_count = list(q)[0]['count']
    assert actual_count == 2158, 'Expected from:Gene, chromosome:22 to:Transcript expected:2158 actual: {}'.format(actual_count)
