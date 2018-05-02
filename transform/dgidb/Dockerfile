FROM golang:1.9.3
ADD *.go /opt/
WORKDIR /opt/
RUN go get github.com/biostream/schemas/go/bmeg
RUN go build compound-id-download.go
RUN go build dgidb-download.go
RUN go build dgidb-transform.go
ENV PATH="/opt/:${PATH}"
VOLUME /out
WORKDIR /out
