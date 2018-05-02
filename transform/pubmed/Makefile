
all: bmeg-code

bmeg-code:
	cd nlp-schema && \
	protoc \
	-I . \
	--python_out=../ \
	nlp.proto
