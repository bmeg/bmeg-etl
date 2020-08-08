
rule HUMAN_9606_idmapping_dat_gz:
	input:
		"source/go/version.txt"
	output:
		"source/go/HUMAN_9606_idmapping.dat.gz"
	shell:		"curl --verbose --progress-bar --ipv4 --connect-timeout 8 --max-time 120 --retry 128 --ftp-ssl --disable-epsv --ftp-pasv ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping.dat.gz --output source/go/HUMAN_9606_idmapping.dat.gz"

rule goa_human_gaf_gz:
	input:
		"source/go/version.txt"
	output:
		"source/go/goa_human.gaf.gz"
	shell:		"wget http://www.geneontology.org/gene-associations/goa_human.gaf.gz -O source/go/goa_human.gaf.gz"

rule obo:
	input:
		"source/go/version.txt"
	output:
		"source/go/go.obo"
	shell:		"wget http://purl.obolibrary.org/obo/go.obo -O source/go/go.obo"

rule go:
	input:
		"transform/go/go_obo2schema.py",
		"source/go/go.obo"
	output:
		"outputs/go/GeneOntologyTerm.Vertex.json.gz",
		"outputs/go/GeneOntologyTerm_ChildTerms_GeneOntologyTerm.Edge.json.gz",
		"outputs/go/GeneOntologyTerm_ParentTerms_GeneOntologyTerm.Edge.json.gz"
	shell:		"python3 transform/go/go_obo2schema.py"

rule gaf2schema:
	input:
		"transform/go/go_gaf2schema.py",
		"source/go/goa_human.gaf.gz",
		"source/go/HUMAN_9606_idmapping.dat.gz"
	output:
		"outputs/go/GeneOntologyTerm_Genes_Gene.Edge.json.gz",
		"outputs/go/Gene_GeneOntologyTerms_GeneOntologyTerm.Edge.json.gz"
	shell:		"python3 transform/go/go_gaf2schema.py"
