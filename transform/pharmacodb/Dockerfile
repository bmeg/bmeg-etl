FROM rocker/r-base:4.3.0

#establishing working directory so we can mount easier
WORKDIR /usr/src

#copying over the two R scripts we'd be using
COPY drug_response_tsv.R ./
COPY download_pharmaco.R ./

#installing library dependencies

#adding curl so R libraries get installed ok
RUN apt-get -y update; apt-get -y install curl
RUN apt-get install -y libcurl4-openssl-dev

RUN R -e "install.packages('BiocManager',dependencies=TRUE)"
RUN R -e "BiocManager::install('PharmacoGx')"
RUN R -e "install.packages('optparse',dependencies=TRUE)"
RUN R -e "install.packages('stringr',dependencies=TRUE)"

