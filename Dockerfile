FROM ubuntu:18.04

RUN apt-get update && \
apt-get install -y python python-pip python-lzo zlib1g-dev unzip

RUN pip install numpy bx-python requests protobuf biostream-schemas pandas xlrd

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
