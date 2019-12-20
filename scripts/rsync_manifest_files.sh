#!/bin/bash

set -e

if [ "$#" -ne 2 ]; then
		printf "Illegal number of parameters.\n\n"
		printf "Usage: rsync_manifest_files.sh FILE_MANIFEST USER@HOST:DEST\n"
		exit 1
fi

file_manifest=$(realpath $1)
dst=$2

copyFile() {
		if [ "$#" -ne 2 ]; then
				printf "Illegal number of parameters.\n\n"
				printf "Usage: copyFile SOURCE USER@HOST:DEST\n"
				exit 1
		fi
		local_file=$1
		remote_loc=$2
		remote=$(echo $remote_loc | tr ":" "\n" | head -n 1)
		dest_root=$(echo $remote_loc | tr ":" "\n" | tail -n 1)
		if [ ! -f ${local_file} ]; then
				echo "file ${local_file} does not exist!"
		else
				fname=$(basename $local_file)
				fdir=$(dirname $(echo $local_file | sed s~outputs/~~))
				echo "COPYING: ${local_file} TO ${remote}:${dest_root}/${fdir}/${fname}"
				echo "$remote"
				ssh $remote "mkdir -p ${dest_root}/${fdir}"
				rsync ${local_file} ${remote}:${dest_root}/${fdir}/${fname}
		fi
}

export -f copyFile

rsync $file_manifest ${dst}/

cat $file_manifest | xargs -L 1 -P 10 -I {} bash -c "copyFile {} $dst"
