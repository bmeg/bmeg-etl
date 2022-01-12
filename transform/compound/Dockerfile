FROM ubuntu:18.04

RUN apt-get update && apt-get install -y python-rdkit librdkit1 rdkit-data python-numpy

ADD chembl_extract.py compound_dist.py /opt/
