# cosmic-extract

Download and extract the cosmic dependencies for g2p



## Build

```
$docker build -t cosmic-extract .
```


## container documentation

```
$ docker inspect cosmic-extract    | jq -r '.[0].Config.Labels'
{
  "maintainer": "walsbr@ohsu.edu",
  "org.label-schema.description": "Download cosmic variants and create tsv(s).",
  "org.label-schema.docker.cmd": "docker run --rm -it -v /tmp:/output -v ~/g2p-aggregator/harvester/CosmicMutantExport.tsv.gz:/g2p-aggregator/harvester/CosmicMutantExport.tsv.gz  cosmic-extract",
  "org.label-schema.name": "g2p cosmic",
  "org.label-schema.usage": "https://github.com/biostream/cosmic-extract/blob/master/README.md",
  "version": "1.0"
}
```

### inputs

* `/g2p-aggregator/harvester/CosmicMutantExport.tsv.gz`
  * you will need to register and download https://grch37-cancer.sanger.ac.uk/cosmic/files?data=/files/grch37/cosmic/v81/CosmicMutantExport.tsv.gz


### outputs

* `/output`
* upload to swift

```
$ swift upload --object-name cosmic etl-development /tmp/cosmic
```

### Example CWL invocation
```
cwltool cosmic-extract.cwl test.json
```

### known consumers

* g2p-transform https://github.com/biostream/g2p-transform
