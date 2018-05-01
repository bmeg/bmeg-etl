

all: bmeg-code ga4gh-code

bmeg-code:
	cd phenotype-schema && \
	protoc \
	-I . -I ../ga4gh-schemas/src/main/proto/ \
	--python_out=../ \
	phenotype.proto

ga4gh-code:
	cd ga4gh-schemas/src/main/proto && \
	protoc \
	-I. \
	--python_out=../../../../ \
	ga4gh/genotype_phenotype.proto 	ga4gh/common.proto ga4gh/bio_metadata.proto
