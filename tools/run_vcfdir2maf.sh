#!/bin/bash

VCFDIR=$1
OUTDIR=$2
TMP=`mktemp -d -p ./`

tar -C $TMP -xvzf `pwd`/source/vep/vep_supporting_files.tar.gz

export TMPDIR=$TMP

for f in $VCFDIR/*.vcf; do
    name=$(echo $(basename $f) | sed 's/.vcf//')
    if [ "$name" != "CCLE_DepMap_18q3_maf_20180718" ]; then
	if [ ! -e source/ccle/mafs/$name ]; then
	    cwltool --outdir $OUTDIR/$name tools/vcf2maf-tools/vcf2maf.cwl --inputVCF $VCFDIR/$name.vcf --refFasta $TMPDIR/home/exacloud/lustre1/SpellmanLab/chiotti/gdan_pipelines/mc3_dev/vcfs/vep/homo_sapiens/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz --vepData $TMPDIR/home/exacloud/lustre1/SpellmanLab/chiotti/gdan_pipelines/mc3_dev/vcfs/vep
	fi
    fi
done

rm -rf $TMP
