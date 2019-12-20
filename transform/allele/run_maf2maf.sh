#!/bin/bash

set -e

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
				mkdir -p $VEPTMP
				tar -C $VEPTMP -xzf source/vep/vep_supporting_files.tar.gz
		else
				echo "skipping untar; found $VEPTMP"
		fi
		export TMPDIR=`mktemp -d -p ./`
		cwltool --outdir $OUTDIR \
						--leave-tmpdir \
						--copy-outputs \
						tools/vcf2maf-tools/maf2maf.cwl \
						--inputMAF $MINIMAL_MAF \
						--outputMAF $(basename $ANNOTATED_MAF) \
						--vepData $VEPTMP/home/exacloud/lustre1/SpellmanLab/chiotti/gdan_pipelines/mc3_dev/vcfs/vep \
						--refFasta $VEPTMP/home/exacloud/lustre1/SpellmanLab/chiotti/gdan_pipelines/mc3_dev/vcfs/vep/homo_sapiens/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz \
						--vepForks 4 \
						--bufferSize 4000
		# rm -rf $TMPDIR
else
		echo "skipping maf2maf conversion; found $ANNOTATED_MAF"
fi
