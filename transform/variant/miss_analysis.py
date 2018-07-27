#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" summarize logs to analyze misses  """

# print total imported
# cat  maf_transform.log | sed 's/2018.* - .* - //' | grep imported  | tail -1 | sed 's/imported //'

# cat  maf_transform.log | sed 's/2018.* - .* - //' | grep stage | python miss_analysis.py  $(cat  maf_transform.log | sed 's/2018.* - .* - //' | grep imported  | tail -1 | sed 's/imported //') |  jq .
# {
#   "myvariantinfo_nofind": 26302,
#   "myvariantinfo_exception": 16,
#   "dbSNP_mismatch": 2926,
#   "Variant_Type": {
#     "ONP": 13,
#     "DEL": 4582,
#     "TNP": 4,
#     "SNP": 21782,
#     "INS": 2863
#   }
# }

import json
import sys

total = int(sys.argv[1])
summary = {'total': total}


for line in sys.stdin:
    line = eval(line)
    summary[line['stage']] = summary.get(line['stage'], 0) + 1
    for k, v in [a.split(':') for a in line['annotations']]:
        summary[k] = summary.get(k, {})
        summary[k][v] = summary[k].get(v, 0) + 1
    if line['reference_bases'] == '-':
        summary['reference_wildcard'] = summary.get('reference_wildcard', 0) + 1
    if line['alternate_bases'] == '-':
        summary['alternate_wildcard'] = summary.get('alternate_wildcard', 0) + 1

summary['missing_snp'] = summary['dbSNP_RS']['novel'] + summary['dbSNP_RS']['.']

summary['report'] = 'misses = {}%; novel = {}%; wildcard_misses ={}%; dbSNP_mismatch ={}%'.format(
    ((summary['myvariantinfo_nofind']+summary['dbSNP_mismatch']) *100)/summary['total'],
    (summary['missing_snp']*100)/summary['total'],
    ((summary['alternate_wildcard'] + summary['reference_wildcard'])*100)/summary['total'],
    (summary['dbSNP_mismatch']*100)/summary['total'],
)

del summary['Feature']
del summary['dbSNP_RS']
del summary['Feature_type']
print(json.dumps(summary))
