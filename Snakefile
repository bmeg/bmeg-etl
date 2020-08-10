
OUTPUTS = []

include: "Snakefile.ccle"
include: "Snakefile.pfam"
#include: "Snakefile.prism"
include: "Snakefile.ensembl"
include: "Snakefile.tcga"
include: "Snakefile.mc3"
#include: "Snakefile.vep"
#include: "Snakefile.dgidb"
#include: "Snakefile.pubmed"
#include: "Snakefile.go"
#include: "Snakefile.pathway_commons"
#include: "Snakefile.gtex"
#include: "Snakefile.gene_enricher"
#include: "Snakefile.g2p"
#include: "Snakefile.msigdb"
include: "Snakefile.gdc"
#include: "Snakefile.mondo"
#include: "Snakefile.publication"
#include: "Snakefile.allele"
include: "Snakefile.pharmacodb"
include: "Snakefile.ctrp"
include: "Snakefile.gdsc"
include: "Snakefile.celllines"
#include: "Snakefile.phenotype"
#include: "Snakefile.compound"

rule outputs_bmeg_manifest:
	input:
		OUTPUTS
	shell:
		"echo generating file manifest..."
