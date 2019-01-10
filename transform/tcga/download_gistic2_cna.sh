#!/bin/bash

./transform/tcga/external_scripts/firehose_get -b -tasks gistic analyses latest

mkdir -p source/tcga/gistic2-firehose/tmp/

mv ./analyses__2016_01_28 source/tcga/gistic2-firehose/tmp/

find source/tcga/gistic2-firehose/tmp/ -type f -name *.CopyNumber_Gistic2.Level_4.*.tar.gz | sed -e 's~\(source/tcga/gistic2-firehose/tmp/\)\(analyses__2016_01_28/.*/.*/\)\(gdac.broadinstitute.org_\)\(.*\)\(-.*.CopyNumber_Gistic2.Level_4.*tar.gz\)~mkdir \1\4; tar -xzf \1\2\3\4\5 -C \1\4~g' | bash

find source/tcga/gistic2-firehose/tmp/ -type f -name all_thresholded.by_genes.txt | sed -e 's~\(source/tcga/gistic2-firehose/\)\(tmp/\)\(.*\)\(/.*/\)\(all_thresholded.by_genes.txt\)~mv \1\2\3\4\5 \1TCGA-\3_\5~' | bash

rm -r source/tcga/gistic2-firehose/tmp/
