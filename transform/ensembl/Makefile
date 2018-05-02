
GOPATH := $(shell pwd)
export GOPATH
PATH := ${PATH}:$(shell pwd)/bin
export PATH

proto:
	cd genome-schema && \
	protoc \
	$(PROTO_INC) \
	--go_out=../src/bmeg \
	genome.proto

download:
	go get github.com/golang/protobuf/protoc-gen-go
