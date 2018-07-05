FROM python:2.7

# build utilities
RUN apt-get update ; apt-get install build-essential -y
RUN apt-get install unzip gzip  -y


# get source & build it
RUN git clone https://github.com/ohsu-comp-bio/g2p-aggregator.git
WORKDIR /g2p-aggregator
RUN git checkout new-harvesters

# install dependencies
WORKDIR /g2p-aggregator/harvester
RUN pip install -r requirements.txt

# output should appear here
VOLUME ["/output"]

COPY ./make-all.sh /g2p-aggregator/harvester



# documentation
LABEL maintainer="walsbr@ohsu.edu" \
 org.label-schema.name="g2p cosmic" \
 org.label-schema.description="Download cosmic variants and create tsv(s)." \
 org.label-schema.usage="https://github.com/biostream/cosmic-extract/blob/master/README.md" \
 version="1.0"

LABEL org.label-schema.docker.cmd="docker run --rm -it -v /tmp:/output -v ~/g2p-aggregator/harvester/CosmicMutantExport.tsv.gz:/g2p-aggregator/harvester/CosmicMutantExport.tsv.gz  cosmic-extract"


# make file uses gzcat
RUN  ln -s /bin/zcat /bin/gzcat
