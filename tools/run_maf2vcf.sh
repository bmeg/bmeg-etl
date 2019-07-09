#!/bin/bash

MAFFILE=$1
OUTDIR=$2

if [ ! -e $OUTDIR ]; then
		echo "starting maf2vcf conversion..."
		TMP=`mktemp -d -p ./`
		tar -C $TMP -xvzf `pwd`/source/vep/vep_supporting_files.tar.gz
		export TMPDIR=$TMP
		cwltool --outdir $OUTDIR tools/vcf2maf-tools/maf2vcf.cwl --perTnVcfs --inputMaf $MAFFILE --refFasta $TMP/home/exacloud/lustre1/SpellmanLab/chiotti/gdan_pipelines/mc3_dev/vcfs/vep/homo_sapiens/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz
		rm -rf $TMP
else
		echo "skipping maf2vcf conversion; outputs exist in $OUTDIR"
fi
