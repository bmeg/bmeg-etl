

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

Note flake8 has a patch for python3.7, (prints a futures warning) I installed it via 
```
pip install -e git+https://gitlab.com/pycqa/flake8@9631dac5#egg=flake8 --src /usr/local/lib/python3.7/site-packages
```



### Functionality

##### Simple mismatches


Given this entry in a MAF file

```
{'stage': 'dbSNP_mismatch', 'genome': 'GRCh37', 'chromosome': '7', 'start': 107763620, 'end': 107763620, 'reference_bases': 'G', 'alternate_bases': 'T', 'annotations': ['Variant_Type:SNP', 'Feature_type:Transcript', 'Feature:ENST00000388781', 'dbSNP_RS:rs199901974'] }
```

A lookup shows that while there is a hit  for that id `rs199901974` the alternate does not match:

```
$ curl -s  http://myvariant.info/v1/variant/rs199901974 | jq .vcf
{
  "alt": "A",
  "position": "107763620",
  "ref": "G"
}
```
A manual lookup of dbSNP directly agrees with myvariantinfo.
![image](https://user-images.githubusercontent.com/47808/43287246-4d45670c-90d9-11e8-9648-d8ad557e1c24.png)

`Since 'G>T' != 'A>G' we categorize this as a "miss"`



#### Inserts & Deletions

`Should we have a better matching algorithm? How will it impact id hash`


From  https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/
```
18 - Match_Norm_Seq_Allele1	Primary data genotype. Matched normal sequencing allele 1. A "-" symbol for a deletion represents a variant. A "-" symbol for an insertion represents wild-type allele. Novel inserted sequence for insertion does not include flanking reference bases (cleared in somatic MAF)


20 - Tumor_Validation_Allele1	Secondary data from orthogonal technology. Tumor genotyping (validation) for allele 1. A "-" symbol for a deletion represents a variant. A "-" symbol for an insertion represents wild-type allele. Novel inserted sequence for insertion does not include flanking reference bases
```

Given this entry in a MAF file

```
{'genome': 'GRCh37', 'chromosome': '6', 'start': 107531642, 'end': 107531642, 'reference_bases': 'C', 'alternate_bases': '-', 'annotations': ['Variant_Type:DEL', 'Feature_type:Transcript', 'Feature:ENST00000369037', 'dbSNP_RS:rs780560462'] ...
```

A lookup shows that while there is a hit  for that id `rs780560462` the alternate does not match:

```
curl -s  http://myvariant.info/v1/variant/rs780560462 | jq .vcf
{
  "alt": "A",
  "position": "107531642",
  "ref": "C"
}
```


Similarly, the same holds true for a reference_base

```
{'genome': 'GRCh37', 'chromosome': '12', 'start': 54757553, 'end': 54757554, 'reference_bases': '-', 'alternate_bases': 'C', 'annotations': ['Variant_Type:INS', 'Feature_type:Transcript', 'Feature:ENST00000551809', 'dbSNP_RS:rs771777487'] ...
```
A lookup shows that while there is a hit  for that id `rs771777487` the reference does not match:

```
$ curl -s 'http://myvariant.info/v1/variant/rs771777487' | jq .vcf
{
  "alt": "AC",
  "position": "54757553",
  "ref": "A"
}
```
