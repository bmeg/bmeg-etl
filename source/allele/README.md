

## `allele.Allele.Vertex.json`

This file is created by:

```
python3 -m transform.allele.transform
```

## Processing steps:

* Reads the output vertex of other transformers, sorts and dedups them while merging allele.annotations.  Output from this stage is retained at `source/allele/sorted_allele_file.json`
  * outputs/g2p/g2p.Allele.Vertex.json.gz
  * outputs/mc3/mc3.Allele.Vertex.json.gz
  * outputs/ccle/ccle.Allele.Vertex.json.gz

* Decorates them with source/myvariantinfo/biothings_current_old_hg19.json.gz. Output from this stage is retained at `source/allele/sqlite.db`
