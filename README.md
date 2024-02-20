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

Get all files needed for build
```
lathe  run prep.plan
```


Run build
```
lathe run build.plan
```

## Visualize Build

Lathe can render the DAG of operations using GraphViz

```
lathe viz build.plan | dot -Tpng -o test.png
```