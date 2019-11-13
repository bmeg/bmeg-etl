#!/bin/bash

if [ "$#" -ne 2 ]; then
    echo "Illegal number of parameters."
		echo "Usage: run_maf2maf.sh <minimal-maf> <output-filename>"
		exit 1
fi

MINIMAL_MAF=$1
ANNOTATED_MAF=$2
OUTDIR=$(dirname $ANNOTATED_MAF)

if [ ! -e $ANNOTATED_MAF ]; then
		echo "running maf2maf..."
		VEPTMP=source/vep/tmp
		if [ ! -e $VEPTMP ]; then
				tar -C $TMP -xzf $PWD/source/vep/vep_supporting_files.tar.gz
		else
				echo "skipping untar; found $VEPTMP"
		fi
		export TMPDIR=`mktemp -d -p ./`
		cwltool --outdir $OUTDIR tools/vcf2maf-tools/maf2maf.cwl \
						--inputMAF $MINIMAL_MAF \
						--outputMAF $(basename $ANNOTATED_MAF) \
						--vepData $VEPTMP/home/exacloud/lustre1/SpellmanLab/chiotti/gdan_pipelines/mc3_dev/vcfs/vep \
						--refFasta $VEPTMP/home/exacloud/lustre1/SpellmanLab/chiotti/gdan_pipelines/mc3_dev/vcfs/vep/homo_sapiens/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz
		rm -rf $TMPDIR
else
		echo "skipping maf2maf conversion; found $ANNOTATED_MAF"
fi
