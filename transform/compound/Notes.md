
# Notes from creating compound translation tables


## Preliminary scan

 - run `./transform/compound/transform.py`
 - examine contents of `missing/compounds.tsv`
 - Prelim scan creation
```
cat missing/compounds.tsv | grep "^gdc" | awk -F "\t" '{print $2}' | ./transform/compound/normalize.py > source/compound/gdc.prelim
cat missing/compounds.tsv | grep "^prism" | awk -F "\t" '{print $2}' | ./transform/compound/normalize.py > source/compound/prism.prelim
cat missing/compounds.tsv | grep "^g2p" | awk -F "\t" '{print $2}' | ./transform/compound/normalize.py > source/compound/g2p.prelim.table
cat missing/compounds.tsv | grep "^pharmacodb" | awk -F "\t" '{print $2}' | ./transform/compound/normalize.py > source/compound/pharmacodb.prelim
cat missing/compounds.tsv | grep "^dgidb" | awk -F "\t" '{print $2}' | ./transform/compound/normalize.py > source/compound/dgidb.prelim
```

Filter for things with hits and create conversion tables
```
cat source/compound/pharmacodb.prelim | awk -F "\t" '{if (length($3)) print $1 "\t" $2}' > reference/compound/pharmacodb.table
cat source/compound/prism.prelim | awk -F "\t" '{if (length($3)) print $1 "\t" $2}' > reference/compound/prism.table
cat source/compound/g2p.prelim | awk -F "\t" '{if (length($3)) print $1 "\t" $2}' > reference/compound/g2p.table
```

For some reason some of the terms are normalized to DrugBank IDs, but then can't be found. Debug
```
cat source/compound/gdc.prelim | awk -F "\t" '{if (length($3)) print $1 "\t" $2}'  | egrep -v "DB[0-9]*" > reference/compound/gdc.table
cat source/compound/dgidb.prelim | awk -F "\t" '{if (length($3)) print $1 "\t" $2}' | egrep -v "DB[0-9]*" > reference/compound/dgidb.table
cat source/compound/prism.prelim | awk -F "\t" '{if (length($3)) print $1 "\t" $2}' | egrep -v "DB[0-9]*" > reference/compound/prism.table
```
