#!/bin/bash

set -e

file=./source/tcga/methylation/IlluminaHumanMethylation450.tar.gz
fname=`basename $file`
ftype=${fname%%.*}
outdir=./source/tcga/methylation
tmp_dir=${outdir}/${ftype}_tmp
tmp=${outdir}/${ftype}.paste.tmp
tmp_out=${outdir}/${ftype}.matrix.tmp
out=${outdir}/${ftype}.tsv

if [ -f "${out}.gz" ]; then
		echo "$(date +%Y-%m-%d\ %H:%M:%S)  matrix ${out}.gz already exists"
		exit 0
fi

if [ ! -d "$tmp_dir" ]; then
		echo "$(date +%Y-%m-%d\ %H:%M:%S)  unarchiving"
		mkdir $tmp_dir
		tar -C $tmp_dir -xf ${file}
fi

echo "$(date +%Y-%m-%d\ %H:%M:%S)  creating matrix"
sample_files=($(ls ${tmp_dir}/*/*.txt))
colnames=("Composite Element REF")
paste_cmd="paste <(cat $tmp_out)"
for i in ${!sample_files[@]}; do
		echo "$(date +%Y-%m-%d\ %H:%M:%S)  processing ${i}/${#sample_files[@]}"
		sample=$(basename ${sample_files[$i]} .txt | tr "." "\t"| cut -f 6)
		colnames+=("$sample")
		if [ $i -eq 0 ]; then 
				cat ${sample_files[$i]} | sed '1,2d' | cut -f 1,2 > $tmp_out
		elif [ $i -eq $((${#sample_files[@]} - 1)) ]; then
				colnames+=("Gene_Symbol" "Chromosome" "Genomic_Coordinate")
				paste_cmd=${paste_cmd}" <(cat ${sample_files[$i]} | sed '1,2d' | cut -f 2,3,4,5)"				
		else
				paste_cmd=${paste_cmd}" <(cat ${sample_files[$i]} | sed '1,2d' | cut -f 2)"
		fi
		
		if [ $(($i % 200)) == 0 ]; then
				eval $paste_cmd > $tmp
				mv $tmp $tmp_out
				paste_cmd="paste <(cat $tmp_out)"
		fi
done

eval $paste_cmd > $tmp
mv $tmp $tmp_out

echo "$(date +%Y-%m-%d\ %H:%M:%S)  adding header"
printf "%s\t" "${colnames[@]}" > ${tmp}_header
printf "\n" >> ${tmp}_header
cat ${tmp}_header $tmp_out > $tmp
mv $tmp $out
echo "$(date +%Y-%m-%d\ %H:%M:%S)  cleaning up"
rm ${tmp}_header
rm ${tmp_out}
rm -rf ${tmp_dir}
echo "$(date +%Y-%m-%d\ %H:%M:%S)  running gzip"
gzip $out
echo "$(date +%Y-%m-%d\ %H:%M:%S)  done"
