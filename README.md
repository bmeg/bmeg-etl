[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

BMEG-ETL
========

BMEG-ETL is a project that defines the transformers and data models for BMEG.


Compile docker image for development build
```
./runner.py docker-build
```

Running pipelines
-----------------

Download source files
```
mkdir build
cd build
swift download biostream -p source
```

Run single pipelines
```
./runner.py run config/ensembl.cwl.yaml build/
```
