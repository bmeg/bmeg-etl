
"""
a set of very basic queries - simply ensure the counts of label->label
"""

import time
from gripql import eq
import json

EXPECTED_COUNTS = [
    {'_from': 'Biosample', 'to': 'Individual', 'expected_count': 57186},
    {'_from': 'Individual', 'to': 'Project', 'expected_count': 32555},
    {'_from': 'Aliquot', 'to': 'Biosample', 'expected_count': 171813},
    {'_from': 'DrugResponse', 'to': 'Biosample', 'expected_count': 169865},
    {'_from': 'Protein', 'to': 'PFAMFamily', 'expected_count': 108729},
    {'_from': 'Protein', 'to': 'Transcript', 'expected_count': 73439},
    {'_from': 'Callset', 'to': 'Aliquot', 'expected_count': 276975},
    # ? {'_from': 'Allele', 'to': 'Callset', 'expected_count': 1},
]


class Stopwatch:

    # Construct self and start it running.
    def __init__(self):
        self.creationTime = time.time()  # Creation time

    # Return the elapsed time since creation of self, in seconds.
    def elapsedTime(self):
        return time.time() - self.creationTime


def count_traversal(_from, to, expected_count, V, expected_time=60):
    """ count traversal template query """
    watch = Stopwatch()
    q = V.where(eq("_label", _from)).out().where(eq("_label", to)).count()
    query_string = json.dumps(q.__dict__, separators=(',',':'))
    actual_count = list(q)[0]['count']
    actual_time = watch.elapsedTime()
    if actual_count != expected_count:
        return 'Expected from:{} to:{} expected:{} actual: {}'.format(
            _from, to, expected_count, actual_count
        )
    if actual_time > expected_time:
        return 'Time exceeded from:{} to:{} expected:{} actual: {}\n    {}'.format(
            _from, to, expected_time, actual_time, query_string
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
