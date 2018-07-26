

ftp://ftp.ensembl.org/pub/grch37/release-90/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz

ftp://ftp.ensembl.org/pub/grch37/release-90/variation/VEP/homo_sapiens_vep_90_GRCh37.tar.gz


## verify source style
```
# pip install -e git+https://gitlab.com/pycqa/flake8@9631dac5#egg=flake8 --src /usr/local/lib/python3.7/site-packages
flake8  --exclude=ga4gh*.py,tests/*.py
# no output is the expected result
```

## testing
```
python -m pytest   --cov=.
platform linux -- Python 3.7.0, pytest-3.6.3, py-1.5.4, pluggy-0.6.0
rootdir: /src, inifile:
plugins: cov-2.5.1
collected 10 items

tests/test_maf_transform.py ..........                                                                                                                                                                  [100%]

----------- coverage: platform linux, python 3.7.0-final-0 -----------
Name                  Stmts   Miss  Cover   Missing
---------------------------------------------------
__init__.py               0      0   100%
allele_harvester.py      41      1    98%   55
maf_transform.py        103      0   100%
---------------------------------------------------
TOTAL                   144      1    99%

```
