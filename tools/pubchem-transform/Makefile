

all: phenotype-code ga4gh-code

phenotype-code:
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
	ga4gh/genotype_phenotype.proto ga4gh/common.proto 
