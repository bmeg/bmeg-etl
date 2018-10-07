
"""
a set of very basic queries - simply ensure the counts of label->label
"""

import time
import logging

EXPECTED_COUNTS = [
    {'_from': 'Biosample', 'to': 'Individual', 'via': 'BiosampleFor', 'expected_count': 73830, 'expected_time': 22},
    {'_from': 'Individual', 'to': 'Project', 'via': 'InProject', 'expected_count': 34353, 'expected_time': 11},
    {'_from': 'Aliquot', 'to': 'Biosample', 'via': 'AliquotFor', 'expected_count': 187411, 'expected_time': 50},
    {'_from': 'DrugResponse', 'to': 'Biosample', 'via': 'DrugResponseIn', 'expected_count': 169865, 'expected_time': 53},
    {'_from': 'Protein', 'to': 'PFAMFamily', 'via': 'PFAMAlignment', 'expected_count': 108729, 'expected_time': 30},
    {'_from': 'Protein', 'to': 'Transcript', 'via': 'ProteinFor', 'expected_count': 73439, 'expected_time': 29},
    {'_from': 'Callset', 'to': 'Aliquot', 'via': 'CallsetFor', 'expected_count': 276975, 'expected_time': 60},
    # # ? {'_from': 'Allele', 'to': 'Callset', 'expected_count': 1},
]


class Stopwatch:

    # Construct self and start it running.
    def __init__(self):
        self.creationTime = time.time()  # Creation time

    # Return the elapsed time since creation of self, in seconds.
    def elapsedTime(self):
        return time.time() - self.creationTime


def count_traversal(_from, to, expected_count, postgres, via=None, expected_time=60):
    """ count traversal template query """

    sql = """select count(f.from) as "count" from edge as f where f.label = '{}' ; """.format(via)

    watch = Stopwatch()
    actual_count = [x for x in postgres.query(sql)][0]['count']
    actual_time = watch.elapsedTime()
    via_msg = via
    if not via_msg:
        via_msg = 'any'
    print('{} {}-{}->{} {}'.format(actual_time, _from, via_msg, to, actual_count))
    if actual_count != expected_count:
        return 'Expected from:{} to:{} expected:{} actual: {}\n    {}'.format(
            _from, to, expected_count, actual_count, sql
        )
    if actual_time > expected_time:
        return 'Time exceeded from:{} to:{} expected:{} actual: {}\n    {}'.format(
            _from, to, expected_time, actual_time, sql
        )


def test_expected_counts(postgres, caplog):
    """ iterate through EXPECTED_COUNTS, assert expected_count """
    caplog.set_level(logging.INFO)
    errors = []
    for traversal in EXPECTED_COUNTS:
        error_msg = count_traversal(**traversal, postgres=postgres)
        if error_msg:
            errors.append(error_msg)
    if len(errors) != 0:
        for e in errors:
            print(e)
        assert False, 'Should have no errors'


# def test_expected_exon_transcript(postgres, caplog):
#     """ subset only for one chromosome """
#     caplog.set_level(logging.INFO)
#     q = (
#         V
#         .where(eq("_label", 'Gene'))
#         .where(eq("chromosome", '22'))
#         .in_().where(eq("_label", 'Transcript'))
#         .count()
#     )
#     actual_count = list(q)[0]['count']
#     assert actual_count == 2158, 'Expected from:Gene, chromosome:22 to:Transcript expected:2158 actual: {}'.format(actual_count)
