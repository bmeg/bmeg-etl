

ftp://ftp.ensembl.org/pub/grch37/release-90/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz

ftp://ftp.ensembl.org/pub/grch37/release-90/variation/VEP/homo_sapiens_vep_90_GRCh37.tar.gz


## verify source style
```
flake8
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

#### Execution

```
python maf_transform.py --maf_file source/mc3.v0.2.8.PUBLIC.maf.gz --prefix biostream/mc3 --verbose --workers 10  2> maf_transform.log
```

#### Logs analysis

Strip log format prefix

```
cat maf_transform.log | sed 's/^.* - .* - //' | jq  '. | your log parser '
```

Harvest results

```
TOTAL=$(cat  maf_transform.log | grep imported  | tail -1 | sed 's/^.*imported //')
sed 's/^.* - .* - //' maf_transform.log  | grep stage | python miss_analysis.py $TOTAL | jq .

{
  "total": 3600963,
  "alt_off_by_one": 5213,
  "ref_off_by_one": 12708,
  "myvariantinfo_nofind": 860324,
  "Variant_Type": {
    "SNP": 809789,
    "DEL": 168164,
    "INS": 44531,
    "TNP": 23,
    "ONP": 151
  },
  "alternate_wildcard": 165272,
  "reference_wildcard": 43638,
  "dbSNP_mismatch": 95899,
  "myvariantinfo_exception": 66435,
  "missing_snp": 760394,
  "report": "misses = 26.554646632025932%; novel = 21.116406916705337%; wildcard_misses =5.80150365332829%; dbSNP_mismatch =2.6631487188288245%"
}

```


where:

* myvariantinfo_nofind - nothing found either by location or dbSNP
  * reference_wildcard, alternate_wildcard  - the nofind had a wild card (ref or alt was '-')
  * dbSNP_mismatch - a second lookup was attempted by snp, but the alleles returned did not match ref, alt
* myvariantinfo_exception - parsing or network error
* missing_snp - novel
* Variant_Type - for all misses


### dedup

```
# reduce to unique ids
$jq -r -c .gid  < mc3.Allele.Vertex.json  | sort | uniq >mc3.Allele.Vertex.unique_ids
$ wc -l mc3.Allele.Vertex.unique_ids
3088088 mc3.Allele.Vertex.unique_ids
```

### compare against original

```
# get the unique id of everything we harvested
jq -r -c .gid  < mc3.Allele.Vertex.json  | sort | uniq >mc3.Allele.Vertex.unique_ids

# run without harvesting to create a reference set
$python maf_transform.py --maf_file source/mc3.v0.2.8.PUBLIC.maf.gz --prefix reference/mc3 --verbose --workers 10  2> maf_transform.log

$ wc -l test.Allele.Vertex.unique_ids
3091454 test.Allele.Vertex.unique_ids
```

### compare against g2p

```
# g2p_analysis.py
{"total":46575,"hits":29488,"keys":{"_id":24749,"_score":24749,"chrom":24749,"clinvar":10389,"hg19":24749,"observed":16746,"snpeff":24749,"vcf":24749,"cadd":23517,"dbsnp":13995,"gnomad_genome":4209,"wellderly":4381,"dbnsfp":11918,"evs":3024,"exac_nontcga":4566,"gnomad_exome":8477,"mutdb":7164,"cosmic":10433,"exac":5257,"grasp":203,"geno2mp":3006,"emv":860,"snpedia":276,"gwassnps":21,"civic":4853,"cgi":2943,"docm":4332}}
```

#### myvariant meta coverage

[snpeff](http://snpeff.sourceforge.net/) 24749
[cadd](http://cadd.gs.washington.edu/home) 23517
[dbsnp](https://www.ncbi.nlm.nih.gov/projects/SNP/) 13995
[dbnsfp](https://sites.google.com/site/jpopgen/dbNSFP) 11918
[cosmic](http://cancer.sanger.ac.uk/cosmic) 10433
[clinvar](https://www.ncbi.nlm.nih.gov/clinvar/) 10389
[mutdb](http://www.mutdb.org/) 7164
[exac](http://exac.broadinstitute.org/) 5257
[civic](https://civic.genome.wustl.edu/home) 4853
[wellderly](https://genomics.scripps.edu/browser/) 4381
[docm](None) 4332
[evs](http://evs.gs.washington.edu/EVS/) 3024
[geno2mp](http://geno2mp.gs.washington.edu) 3006
[cgi](https://www.cancergenomeinterpreter.org/home) 2943
[emv](http://www.egl-eurofins.com/emvclass/emvclass.php) 860
[snpedia](https://www.snpedia.com/) 276
[grasp](https://grasp.nhlbi.nih.gov/Updates.aspx) 203
[gwassnps](http://www.ebi.ac.uk/gwas/) 21

#### G2P dbnsfp.Vest3 rankscore

* A 0.8763147826086949
* B 0.6735934167140423
* C 0.6866358542713499
* D 0.8241483648428666
* None 0.784781306053812


#### G2P cadd.rawscore
A 4.8880969561332135
B 1.29712009661252
C 1.3459371859442335
D 2.8110570133636186
None 3.1775750201342317
