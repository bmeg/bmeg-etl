FROM ubuntu:18.04

RUN apt-get update && \
apt-get install -y python3 python3-pip python-lzo zlib1g-dev unzip git

RUN pip install numpy bx-python requests protobuf pandas xlrd
RUN pip install "git+https://github.com/bmeg/schemas.git#subdirectory=python"

# GDC Transform
COPY transform/gdc/*.py /opt/

# Ensembl Gene / Transcript / Exon
COPY transform/ensembl/*.py /opt/

WORKDIR /opt
