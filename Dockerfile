FROM ubuntu:18.04

RUN apt-get update && \
apt-get install -y python python-pip python-lzo zlib1g-dev unzip golang-go git

RUN pip install numpy bx-python requests protobuf pandas xlrd
RUN pip install "git+https://github.com/bmeg/schemas.git#subdirectory=python"

# GDC Transform
COPY transform/gdc/*.py /opt/
COPY transform/gdc/tcga_pubchem.map /opt/

# CCLE Transform
COPY transform/ccle/*.py /opt/
COPY transform/ccle/ccle_pubchem.txt /opt/

# CTDD transform
COPY transform/ctdd/*.py /opt/
COPY transform/ctdd/ctdd_pubchem.table /opt/

# Variant transform
COPY transform/variant/ga4gh-variant.py /opt/

# GO Transform
COPY transform/go/*.py /opt/

# GDSC transform
COPY transform/gdsc/*.py /opt/
COPY transform/gdsc/gdsc_pubchem.table /opt/

COPY transform/ensembl/run.go /opt/

COPY transform/pathway-commons/sif_convert.py /opt/
COPY transform/pubmed/pubmed.py /opt/

ENV GOPATH /opt
RUN go get github.com/blachlylab/gff3
RUN go get github.com/golang/protobuf/jsonpb
RUN go get github.com/golang/protobuf/proto
RUN go get github.com/bmeg/schemas/go/bmeg
