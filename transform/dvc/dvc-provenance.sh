echo  $(dvc pipeline show outputs.bmeg_manifest.dvc | grep source ) $(dvc pipeline show outputs.bmeg_manifest.dvc | grep outputs) | sed 's/ /\n/g' | python transform/dvc/dvc2cmd.py
