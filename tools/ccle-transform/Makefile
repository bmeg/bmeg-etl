

all: phenotype-code rna-code ga4gh-code

phenotype-code:
	cd phenotype-schema && \
	protoc \
	-I . -I ../ga4gh-schemas/src/main/proto/ \
	--python_out=../ \
	phenotype.proto

rna-code:
	cd rna-schema && \
	protoc \
	-I . -I ../ga4gh-schemas/src/main/proto/ \
	--python_out=../ \
	rna.proto

ga4gh-code:
	cd ga4gh-schemas/src/main/proto && \
	protoc \
	-I. \
	--python_out=../../../../ \
	ga4gh/genotype_phenotype.proto 	ga4gh/common.proto ga4gh/bio_metadata.proto ga4gh/metadata.proto \
	ga4gh/variants.proto 	ga4gh/allele_annotations.proto