

all: bmeg-code 

bmeg-code:
	cd pfam-schema && \
	protoc \
	-I . \
	--python_out=../ \
	pfam.proto

