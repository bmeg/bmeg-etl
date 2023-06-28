[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Build Status](https://travis-ci.org/bmeg/bmeg-etl.svg?branch=master)](https://travis-ci.org/bmeg/bmeg-etl)

BMEG-ETL
========

BMEG-ETL is a project that defines the transformers for BMEG. The schema is described here: https://github.com/bmeg/bmeg-dictionary


## Toolchain

### Install GO
Install go following instructions at https://go.dev/dl/

### Install tools
```
go install github.com/bmeg/sifter@main
go install github.com/bmeg/lathe@main
```
### Install other dependencies
```
conda install -c conda-forge rdkit
conda install -c conda-forge pandas
conda install -c conda-forge pyarrow
conda install -c conda-forge polars
```

### Update path

```
export PATH=$PATH:$HOME/go/bin
```

## Building

Build the scripts for generating graph files.
```
lathe plan-graph transform -C graph-build -o graph
```

Build Snakefile
```
lathe plan transform graph-build -C .
```

Run build
```
snakemake -j 4 --resources mem_mb=30000
```

## Loading

### GRIP Schema Generation (note: section needs testing/updates)
------
Use the `generate-schema` executable in the `scripts/` directory to create a schema for the graph. 
This program will dump the graph schema in JSON or YAML format to stdout.

You can rebuild `generate-schema` by running `GOARCH=amd64 GOOS=linux go build -o scripts/generate-schema-linux scripts/generate-schema.go`. 
See https://golang.org/cmd/go/#hdr-Generate_Go_files_by_processing_source for more details or if you need to compile the 
program for another system.

```
$ ./scripts/generate-schema-linux --help
Usage:
  generate-schema [flags]

  Flags:
    -n, --graph-name string       name of the graph
    -h, --help                help for generate-schema
    -m, --manifest string     file manifest listing vertex and edge files
    -s, --sample int          number of elements to sample from each file (default 100)
    -v, --verbose             turn on verbose logging
    -w, --workers int         number of workers to use to read the files and determine the schema (default 10)
    --yaml                    output schema in YAML rather than JSON format
```
