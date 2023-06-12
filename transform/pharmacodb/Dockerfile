FROM rocker/r-base:4.2.1

#establishing working directory so we can mount easier
WORKDIR /src

#copying over the two R scripts we'd be using
COPY ./drug_response_tsv.R /
COPY ./download_pharmaco.R /

#and making sure we are allowed to run them in docker
RUN chmod +x /drug_response_tsv.R /
RUN chmod +x /drug_response_tsv.R /

#installing library dependencies
RUN R -e "install.packages('BiocManager',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('optparse',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('stringr',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "BiocManager::install('PharmacoGx')"

#Preparing code and default arguments for running
ENTRYPOINT ["Rscript"]
CMD ["drug_response_tsv.R"]


