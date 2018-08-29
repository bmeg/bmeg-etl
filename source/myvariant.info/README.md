This directory should contain [
`biothings_current_old_hg19.json.gz`,
`metadata.fields.json`,
`metadata.json`
].

# Download instructions


## create biothings docker container
```
wget http://biothings-containers.s3-website-us-west-2.amazonaws.com/old_myvariant/old_myvariant.docker
docker image load <  old_myvariant.docker
docker run --name old_myvariant -p 19080:80 -p 19200:9200 -p 19022:7022 -p 19090:7080 -d old_myvariant
```
## use biothings container to download mygene.info

see http://docs.biothings.io/en/latest/doc/standalone.html#updating-data-using-biothings-hub

verify results

```
$curl http://localhost:19200/_cat/indices
green  open biothings_current 1 0 14903 0 10.3mb 10.3mb
$curl http://localhost:19080/metadata
{
  "app_revision": ...

```

## use python utility to extract elastic instance as plain old json

`requirements.txt`

```
elasticsearch-dsl>=5.0.0,<6.0.0
```


`get_index.py`

```
#!/usr/bin/python
from elasticsearch import Elasticsearch
from elasticsearch_dsl import Search
import json
import os
import argparse


# Use DSL to query genomic location, subset of fields,
def _to_stdout(index='association'):
    # get connection info from env
    HOST = [os.environ.get('ES', 'localhost:9200')]
    client = Elasticsearch(HOST)
    # validate connection
    assert(client.info()['version'])
    s = Search(using=client, index=index).params(size=1000)
    for hit in s.scan():
        print json.dumps(hit.to_dict(), separators=(',', ':'))



if __name__ == '__main__':
    argparser = argparse.ArgumentParser()
    argparser.add_argument('--index',
                           help='index to write to stdout',
                           )
    args = argparser.parse_args()
    _to_stdout(args.index)
```

Create json
```
export ES=localhost:19200 python get_index.py --index biothings_current > mygene.info.json
curl http://localhost:19080/metadata > mygene.info.metadata.json
```
