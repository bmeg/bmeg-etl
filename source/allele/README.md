
## `minimal_allele.maf`, `minimal_allele.vep.vcf` & `annotated_allele.maf`

These files are created by:

```
python3 transform/allele/harmonize_alleles.py
```

## Processing steps:

* Reads the output vertex of other transformers, sorts and dedups them and then annotates them by running `maf2maf.pl` from https://github.com/mskcc/vcf2maf.
  * outputs/g2p/g2p.Allele.Vertex.json.gz
  * outputs/mc3/mc3.Allele.Vertex.json.gz
  * outputs/ccle/maf.Allele.Vertex.json.gz
  * outputs/gdsc/caveman.Allele.Vertex.json.gz
  * outputs/gdsc/pindel.Allele.Vertex.json.gz
