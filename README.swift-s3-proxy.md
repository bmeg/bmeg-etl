# Usage: s3->swift


## get access token

* after logging into openstack, get access tokens

```
$ openstack ec2 credentials create

+------------+----------------------------------+
| Field      | Value                            |
+------------+----------------------------------+
| access     | AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA |
| project_id | xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx |
| secret     | SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS |
| trust_id   | None                             |
| user_id    | xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx |
+------------+----------------------------------+
```

## aws cli config

* create a new profile in `~/.aws/config`

```
[s3Proxy]
output = json
```

* create a new profile in `~/.aws/credentials`

```
[s3Proxy]
aws_access_key_id           = AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
aws_secret_access_key       = SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
region                      = RegionOne
s3 =
    signature_version = s3
```

note signature_version is critical, swift does not support v3-v4


## verify aws cli config

note `--endpoint` mandatory

```
$ export AWS_PROFILE=s3Proxy
$ aws s3 ls s3://bmeg-dvc    --endpoint=http://10.96.11.20:8080
```


## verify dvc

* push the data to the repository

```
$ ls -lh source/gene_enricher/non_alt_loci_set.json
-rw-r--r--  3 walsbr  OHSUM01\Domain Users    29M Sep 30 21:44 source/gene_enricher/non_alt_loci_set.json

$ dvc remote list
bmeg	s3://bmeg-dvc

$ dvc push -r bmeg
Preparing to push data to s3://bmeg-dvc
[##############################] 100% Collecting information
[##############################] 100% source/gene_enricher/non_alt_loci_set.json

```

* list the repo file

```
$ aws s3 ls s3://bmeg-dvc   --recursive --human-readable --summarize   --endpoint=http://10.96.11.20:8080

2018-10-08 12:41:37   29.4 MiB 06/373d1779ee7f371ced423367ea3144

Total Objects: 1
   Total Size: 29.4 MiB
```

* which corresponds to the dvc local cache

```
$ ls -l .dvc/cache/06/373d1779ee7f371ced423367ea3144
-rw-r--r--  3 walsbr  OHSUM01\Domain Users  30854925 Sep 30 21:44 .dvc/cache/06/373d1779ee7f371ced423367ea3144

```
