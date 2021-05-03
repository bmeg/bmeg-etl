
OUTPUTS = []

COMPOUND_FILES = {}
NORMALIZE_FILES = []
ALLELE_FILES = []

include: "Snakefile.ccle"
include: "Snakefile.pfam"
include: "Snakefile.ensembl"
include: "Snakefile.tcga"
include: "Snakefile.mc3"
include: "Snakefile.dgidb"
include: "Snakefile.go"
include: "Snakefile.pathway_commons"
include: "Snakefile.gtex"
include: "Snakefile.g2p"
include: "Snakefile.msigdb"
include: "Snakefile.gdc"
include: "Snakefile.mondo"
# include: "Snakefile.drugbank"
include: "Snakefile.chembl"

# literature
include: "Snakefile.publication"
include: "Snakefile.pubmed"

include: "Snakefile.allele"
# cellline testing
include: "Snakefile.pharmacodb"
include: "Snakefile.prism"
include: "Snakefile.gdsc"
include: "Snakefile.celllines"
include: "Snakefile.phenotype"
include: "Snakefile.compound"

include: "Snakefile.normalize"

rule outputs_bmeg_manifest:
	input:
		OUTPUTS
	output:
		"bmeg_file_manifest.txt"
	shell:
		"./gen_manifest.sh {input}"
