

all: rna-code clinical-code cna-code ga4gh-code

rna-code:
	cd rna-schema && \
	protoc \
	-I . -I ../ga4gh-schemas/src/main/proto/ \
	--python_out=../ \
	rna.proto

cna-code:
	cd genome-schema && \
	protoc \
	-I . -I ../ga4gh-schemas/src/main/proto/ \
	--python_out=../ \
	cna.proto

clinical-code:
	cd clinical-schema && \
	protoc \
	-I . -I ../ga4gh-schemas/src/main/proto/ \
	--python_out=../ \
	clinical.proto


ga4gh-code:
	cd ga4gh-schemas/src/main/proto && \
	protoc \
	-I. \
	--python_out=../../../../ \
	ga4gh/genotype_phenotype.proto 	ga4gh/common.proto ga4gh/bio_metadata.proto ga4gh/metadata.proto
