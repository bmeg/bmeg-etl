FROM python:2.7

RUN apt-get update && apt-get install -y autoconf automake libtool curl make g++ unzip
WORKDIR /opt
RUN git clone https://github.com/google/protobuf.git && cd protobuf && ./autogen.sh && ./configure && make && make install
RUN ldconfig
RUN pip install protobuf
RUN pip install git+https://github.com/bmeg/bmeg-schemas.git@packaging

