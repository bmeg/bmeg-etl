
"""
a set of very basic queries - simply ensure the counts of label->label
"""

import time
from gripql import eq
import json
import logging

EXPECTED_COUNTS = [
    {'_from': 'Biosample', 'to': 'Individual', 'via': 'BiosampleFor', 'expected_count': 73830, 'expected_time': 22},
    {'_from': 'Individual', 'to': 'Project', 'via': 'InProject', 'expected_count': 34353, 'expected_time': 11},
    {'_from': 'Aliquot', 'to': 'Biosample', 'via': 'AliquotFor', 'expected_count': 187411, 'expected_time': 50},
    {'_from': 'DrugResponse', 'to': 'Biosample', 'via': 'DrugResponseIn', 'expected_count': 169865, 'expected_time': 53},
    {'_from': 'Protein', 'to': 'PFAMFamily', 'via': 'PFAMAlignment', 'expected_count': 108729, 'expected_time': 30},
    {'_from': 'Protein', 'to': 'Transcript', 'via': 'ProteinFor', 'expected_count': 73439, 'expected_time': 29},
    {'_from': 'Callset', 'to': 'Aliquot', 'via': None, 'expected_count': 276975, 'expected_time': 60},
    # # ? {'_from': 'Allele', 'to': 'Callset', 'expected_count': 1},
]

"""
# no edgename
Time Path count
20.152465343475342 Biosample-any->Individual 73830
8.606142044067383 Individual-any->Project 34353
47.062398195266724 Aliquot-any->Biosample 187411
75.5369393825531 DrugResponse-any->Biosample 169865
39.80154490470886 Protein-any->PFAMFamily 108729
40.92285227775574 Protein-any->Transcript 73439
56.05605459213257 Callset-any->Aliquot 276975

# with edgename
~ same, + better, - worse
~ 19.639927625656128 Biosample-BiosampleFor->Individual 73830
- 9.525224924087524 Individual-InProject->Project 34353
~ 48.195698976516724 Aliquot-AliquotFor->Biosample 187411
+ 51.33423638343811 DrugResponse-DrugResponseIn->Biosample 169865
+ 28.148258924484253 Protein-PFAMAlignment->PFAMFamily 108729
+ 27.663743495941162 Protein-ProteinFor->Transcript 73439
- 59.28540301322937 Callset-CallsetFor->Aliquot 276975
"""


class Stopwatch:

    # Construct self and start it running.
    def __init__(self):
        self.creationTime = time.time()  # Creation time

    # Return the elapsed time since creation of self, in seconds.
    def elapsedTime(self):
        return time.time() - self.creationTime


def count_traversal(_from, to, expected_count, V, via=None, expected_time=60):
    """ count traversal template query """
    watch = Stopwatch()
    if via:
        q = V.where(eq("_label", _from)).out(via).where(eq("_label", to)).count()
    else:
        q = V.where(eq("_label", _from)).out().where(eq("_label", to)).count()
    query_string = json.dumps(q.__dict__, separators=(',', ':'))
    actual_count = list(q)[0]['count']
    actual_time = watch.elapsedTime()
    via_msg = via
    if not via_msg:
        via_msg = 'any'
    print('{} {}-{}->{} {}'.format(actual_time, _from, via_msg, to, actual_count))
    if actual_count != expected_count:
        return 'Expected from:{} to:{} expected:{} actual: {}\n    {}'.format(
            _from, to, expected_count, actual_count, query_string
        )
    if actual_time > expected_time:
        return 'Time exceeded from:{} to:{} expected:{} actual: {}\n    {}'.format(
            _from, to, expected_time, actual_time, query_string
        )


def test_expected_counts(V, caplog):
    """ iterate through EXPECTED_COUNTS, assert expected_count """
    caplog.set_level(logging.INFO)
    errors = []
    for traversal in EXPECTED_COUNTS:
        error_msg = count_traversal(**traversal, V=V)
        if error_msg:
            errors.append(error_msg)
    if len(errors) != 0:
        print(errors)
        assert False, 'Should have no errors'


def test_expected_exon_transcript(V, caplog):
    """ subset only for one chromosome """
    caplog.set_level(logging.INFO)
    q = (
        V
        .where(eq("_label", 'Gene'))
        .where(eq("chromosome", '22'))
        .in_().where(eq("_label", 'Transcript'))
        .count()
    )
    actual_count = list(q)[0]['count']
    assert actual_count == 2158, 'Expected from:Gene, chromosome:22 to:Transcript expected:2158 actual: {}'.format(actual_count)
