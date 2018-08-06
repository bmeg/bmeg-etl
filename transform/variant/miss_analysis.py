#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" summarize logs to analyze misses  """

# print total imported
# see README.md for example how to call

import json
import sys
import difflib

total = int(sys.argv[1])
summary = {'total': total, 'alt_off_by_one': 0, 'ref_off_by_one': 0}


for line in sys.stdin:
    try:
        line = eval(line)
    except Exception as e:
        print(line, file=sys.stderr)
        print(e, file=sys.stderr)
        continue

    summary[line['stage']] = summary.get(line['stage'], 0) + 1
    for k, v in [a.split(':') for a in line['annotations']]:
        summary[k] = summary.get(k, {})
        summary[k][v] = summary[k].get(v, 0) + 1
    if line['reference_bases'] == '-':
        summary['reference_wildcard'] = \
            summary.get('reference_wildcard', 0) + 1
    if line['alternate_bases'] == '-':
        summary['alternate_wildcard'] = \
            summary.get('alternate_wildcard', 0) + 1
    if line['stage'] == 'dbSNP_mismatch':
        # print(line['hit'].keys())
        # print(line['hit']['vcf'])
        ref = line['hit']['vcf']['ref']
        alt = line['hit']['vcf']['alt']
        if len(ref) > 1:
            ref_diff = sum([{' ': 0, '-': 1,  '+': 1}[d[0]] for d in difflib.ndiff(line['reference_bases'], ref)])  # noqa
            if ref_diff == 1:
                summary['ref_off_by_one'] += 1
        if len(alt) > 1:
            alt_diff = sum([{' ': 0, '-': 1, '+': 1}[d[0]] for d in difflib.ndiff(line['alternate_bases'], alt)])  # noqa
            if alt_diff == 1:
                summary['alt_off_by_one'] += 1

summary['missing_snp'] = summary['dbSNP_RS'].get('novel', 0) \
                         + summary['dbSNP_RS'].get('.', 0)

report = 'misses = {}%; novel = {}%; wildcard_misses ={}%; dbSNP_mismatch ={}%'
summary['report'] = report.format(
    ((summary['myvariantinfo_nofind']+summary['dbSNP_mismatch']) * 100)/summary['total'],  # noqa
    (summary['missing_snp']*100)/summary['total'],
    ((summary['alternate_wildcard'] + summary['reference_wildcard']) * 100)/summary['total'],  # noqa
    (summary['dbSNP_mismatch']*100)/summary['total'],
)

del summary['Feature']
del summary['dbSNP_RS']
del summary['Feature_type']
print(json.dumps(summary))
