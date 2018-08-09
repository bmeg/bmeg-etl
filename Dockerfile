FROM ubuntu:18.04

RUN apt-get update && \
apt-get install -y python3.7 python3-pip python-lzo zlib1g-dev unzip git liblzo2-dev

RUN pip3 install numpy bx-python requests protobuf pandas xlrd
#RUN pip3 install "git+https://github.com/bmeg/schemas.git#subdirectory=python"

# GDC Transform
COPY transform/gdc/*.py /opt/

# Ensembl Gene / Transcript / Exon
COPY transform/ensembl/*.py /opt/

COPY setup.py README.md /build/
COPY src /build/src

RUN cd /build && python3.7 setup.py build && python3.7 setup.py install

WORKDIR /opt
