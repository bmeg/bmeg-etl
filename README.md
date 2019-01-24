[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Build Status](https://travis-ci.org/bmeg/bmeg-etl.svg?branch=master)](https://travis-ci.org/bmeg/bmeg-etl)

BMEG-ETL
========

BMEG-ETL is a project that defines the transformers and data models for BMEG.


DVC
-----

See:

* [Configure DVC](https://dvc.org/doc/get-started/configure)
* [AWS CLI with Minio Server](https://docs.minio.io/docs/aws-cli-with-minio.html)


Note: All dvc files are maintained in the zip file dvc-graph.zip.  This is due to:

 * each transformer is it's own command, executed with CWD at the project root
 * `dvc run` captures state in a separate \*.dvc file
 * therefore, there are quite a few files at the root level of the project


Setup
-----

```
# until DNS setup, make sure minio.compbio.ohsu.edu is known
sudo sh -c "echo 10.50.50.118 minio.compbio.ohsu.edu >> /etc/hosts"

# see minio install for credentials
# cat /mnt/minio/.minio.sys/config/config.json  | jq .credential

# install and configure aws
$ sudo apt  install awscli
$ pip install awscli
$ aws configure
AWS Access Key ID [None]: KKKKKKKKKKK
AWS Secret Access Key [None]: SSSSSSS
Default region name [None]: us-east-1
Default output format [None]:
# test
aws --endpoint-url https://minio.compbio.ohsu.edu s3 ls
2018-10-29 22:28:14 bmeg


# Setting up minio client
# linux
# install in home directory
cd
wget https://dl.minio.io/client/mc/release/linux-amd64/mc
chmod +x mc
alias mc=~/mc
mc version
# update your config ... vi ~/.mc/config.json
# test
$ mc ls -r  bmeg/bmeg/dvc | head -5
[2018-10-31 00:19:19 UTC] 1.2MiB 07/b930da26e4a06dcf8c9a0faff57be1
[2018-10-31 17:53:28 UTC] 2.6GiB 08/ec6eb40ad76b48210aa0e939ae7aa1
[2018-10-31 00:03:17 UTC] 222MiB 10/a14a8b317a34784e5b3e62c3fa387a
[2018-10-30 23:59:22 UTC]  38MiB 15/525092ad0e0598b95b874c9660bf6c
[2018-10-30 20:23:28 UTC] 809MiB 1c/4711bb30e668d5f387e1819bae99ef

# macOS see brew install

# dvc already installed and initialized
# add our remote
dvc remote add -d minio s3://bmeg/dvc
dvc remote modify minio endpointurl https://minio.compbio.ohsu.edu
```

Example
----------

```
# retrieve data from foreign source
dvc run \
  -o source/gene_enricher/hgnc_complete_set.json \
  --file source.gene_enricher.hgnc_complete_set.json.dvc \
  --yes \
  curl --verbose --progress-bar --ipv4 --connect-timeout 8 --max-time 120 --retry 128 --ftp-ssl --disable-epsv --ftp-pasv ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/json/hgnc_complete_set.json --output source/gene_enricher/hgnc_complete_set.json

# commit DVC's record to github
git add .gitignore source.gene_enricher.hgnc_complete_set.json.dvc
git commit -m "add hgnc_complete_set to DVC"

# commit the data to DVC's remote
$ dvc push
Preparing to push data to s3://bmeg/dvc
[##############################] 100% Collecting information
[##############################] 100% source/gene_enricher/hgnc_complete_set.json

# view the remote repository
$ mc ls -r  bmeg/bmeg/dvc
[2018-10-30 18:32:15 UTC]  29MiB f4/843dade6933b9879654417c6d93c1b

# note that the file name ~ the md5 hash
$ cat source.gene_enricher.hgnc_complete_set.json.dvc
cmd: curl --verbose --progress-bar --ipv4 --connect-timeout 8 --max-time 120 --retry
  128 --ftp-ssl --disable-epsv --ftp-pasv ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/json/hgnc_complete_set.json
  --output source/gene_enricher/hgnc_complete_set.json
md5: 361717476a2235fe752e4c52c774caa1
outs:
- cache: true
  md5: 8300e43e6513e8e0a696a952ad28b1f5
  path: source/gene_enricher/hgnc_complete_set.json
```

Provenance
------

To generate [bmeg_file_manifest.txt](scripts/bmeg_file_manifest.txt) run `bash scripts/generate_bmeg_file_manifest.sh`.

The dvc files were created using the commands in [dvc_commands.txt](dvc_commands.txt)
To recreate this file, run `bash scripts/generate_dvc_commands.sh`.
