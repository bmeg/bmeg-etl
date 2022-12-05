


rule all:
	input:
		"source/pathway_commons/PathwayCommons12.inoh.complex",
		"graph/ensembl.genes-graph.vertex.json.gz",
		"output/mc3/mc3.transform.raw.json.gz",
		"graph/pharmacodb_profiles.cellSample-graph.edge.json.gz",
		"source/pathway_commons/PathwayCommons12.msigdb.extSIF",
		"source/pathway_commons/PathwayCommons12.reactome.complex",
		"source/pathway_commons/PathwayCommons12.innatedb.complex",
		"source/pathway_commons/PathwayCommons12.dip.complex",
		"output/go/go.edgeOpen.edges.json.gz",
		"source/pathway_commons/PathwayCommons12.hprd.extSIF",
		"graph/pubmed.transform-graph.vertex.json.gz",
		"graph/pubmed.transform-graph.edge.json.gz",
		"source/pathway_commons/PathwayCommons12.mirtarbase.complex",
		"output/ensembl/uniprot.transform.protein.json.gz",
		"source/gdsc/model_list_20191104.csv",
		"source/pathway_commons/PathwayCommons12.kegg.complex",
		"output/pharmacodb/sources.tsv.gz",
		"graph/pfam_family.transform-graph.vertex.json.gz",
		"source/pathway_commons/PathwayCommons12.netpath.extSIF",
		"graph/ensembl.transcripts-graph.vertex.json.gz",
		"source/pathway_commons/PathwayCommons12.netpath.complex",
		"source/pathway_commons/PathwayCommons12.drugbank.complex",
		"output/pharmacodb/source_statistics.tsv.gz",
		"graph/pharmacodb_profiles.cellProject-graph.vertex.json.gz",
		"output/pharmacodb/source_tissue_names.tsv.gz",
		"graph/ensembl.genes-graph.edge.json.gz",
		"graph/gdc.caseObject-graph.edge.json.gz",
		"graph/mc3.callset-graph.vertex.json.gz",
		"graph/pharmacodb_profiles.cellCase-graph.edge.json.gz",
		"source/pathway_commons/PathwayCommons12.ctd.extSIF",
		"source/pathway_commons/PathwayCommons12.panther.complex",
		"source/pathway_commons/PathwayCommons12.pathbank.complex",
		"graph/gdc.sampleObject-graph.edge.json.gz",
		"source/ensembl/uniprotId2ensemblGene.alt.tsv",
		"output/gdc/rnaseq.rna.debug.json.gz",
		"source/pathway_commons/PathwayCommons12.humancyc.extSIF",
		"graph/ensembl.exons-graph.edge.json.gz",
		"graph/gdc.aliquotObject-graph.vertex.json.gz",
		"source/pathway_commons/PathwayCommons12.corum.extSIF",
		"source/pathway_commons/PathwayCommons12.reconx.complex",
		"source/pathway_commons/PathwayCommons12.pathbank.extSIF",
		"graph/ensembl.exons-graph.vertex.json.gz",
		"graph/msigdb.object-graph.vertex.json.gz",
		"graph/msigdb.object-graph.edge.json.gz",
		"source/gdc/files.json",
		"output/mondo/mondo.extractNodes.nodes.json.gz",
		"source/pathway_commons/PathwayCommons12.corum.complex",
		"graph/pharmacodb_profiles.cellSample-graph.vertex.json.gz",
		"output/msigdb/msigdb.transform.raw.json.gz",
		"source/pathway_commons/PathwayCommons12.inoh.extSIF",
		"source/pathway_commons/PathwayCommons12.humancyc.complex",
		"source/pathway_commons/PathwayCommons12.reactome.extSIF",
		"output/pharmacodb/dataset_cells.tsv.gz",
		"graph/mc3.callset-graph.edge.json.gz",
		"source/pathway_commons/PathwayCommons12.psp.extSIF",
		"output/pharmacodb/source_cell_names.tsv.gz",
		"graph/pharmacodb_profiles.cellAliquot-graph.edge.json.gz",
		"output/mc3/mc3.variant.allele.json.gz",
		"output/pharmacodb/tissues.tsv.gz",
		"output/pharmacodb/source_drug_names.tsv.gz",
		"graph/chembDrugMechanismExtract.build-graph.edge.json.gz",
		"output/mondo/mondo.extractEdges.edges.json.gz",
		"source/pathway_commons/PathwayCommons12.innatedb.extSIF",
		"output/pharmacodb/cellosaurus.tsv.gz",
		"output/pharmacodb/cell_tissues.tsv.gz",
		"output/go/goReduce.edgeReduce.edge_reduce.json.gz",
		"source/pathway_commons/PathwayCommons12.drugbank.extSIF",
		"source/pathway_commons/PathwayCommons12.dip.extSIF",
		"graph/gdc.caseObject-graph.vertex.json.gz",
		"graph/gdc.sampleObject-graph.vertex.json.gz",
		"source/pathway_commons/PathwayCommons12.bind.extSIF",
		"source/pathway_commons/PathwayCommons12.kegg.extSIF",
		"source/pathway_commons/PathwayCommons12.hprd.complex",
		"graph/go.nodeOpen-graph.vertex.json.gz",
		"graph/pharmacodb_profiles.drObject-graph.edge.json.gz",
		"output/go/go_gaf.dump.gaf.json.gz",
		"graph/pharmacodb_profiles.cellProject-graph.edge.json.gz",
		"source/pathway_commons/PathwayCommons12.msigdb.complex",
		"graph/gdc.aliquotObject-graph.edge.json.gz",
		"graph/go.nodeOpen-graph.edge.json.gz",
		"source/pathway_commons/PathwayCommons12.biogrid.extSIF",
		"source/pathway_commons/PathwayCommons12.mirtarbase.extSIF",
		"source/pathway_commons/PathwayCommons12.ctd.complex",
		"source/pathway_commons/PathwayCommons12.biogrid.complex",
		"source/pathway_commons/PathwayCommons12.panther.extSIF",
		"output/pfam/pfam_family.transform.debug.json.gz",
		"output/pharmacodb/drug_annots.tsv.gz",
		"output/pharmacodb/pharmacodb_profiles.cellDistinct.checkpoint.json.gz",
		"graph/pharmacodb_profiles.cellCase-graph.vertex.json.gz",
		"graph/pharmacodb_profiles.cellAliquot-graph.vertex.json.gz",
		"source/pathway_commons/PathwayCommons12.pid.extSIF",
		"graph/ensembl.transcripts-graph.edge.json.gz",
		"graph/pharmacodb_profiles.drObject-graph.vertex.json.gz",
		"source/ccle/sample_info.csv",
		"output/go/go.nodeOpen.before.json.gz",
		"source/pathway_commons/PathwayCommons12.reconx.extSIF",
		"graph/chembDrugMechanismExtract.build-graph.vertex.json.gz",
		"graph/pfam_family.transform-graph.edge.json.gz",
		"output/msigdb/msigdb.transform.record.json.gz",
		"source/pathway_commons/PathwayCommons12.bind.complex",
		"source/pathway_commons/PathwayCommons12.pid.complex",
		"source/pathway_commons/PathwayCommons12.psp.complex"

rule cell_line_names:
	input:
		"transform/source/pharmacodb/cellosaurus.tsv.gz"
	output:
		"source/ccle/cellline_id_lookup.tsv"
	shell:
		"sifter run transform/ccle/cell_lines.yaml"

rule download:
	output:
		"source/ccle/sample_info.csv"
	shell:
		"cd transform/ccle && curl -L -o ../../source/ccle/sample_info.csv https://ndownloader.figshare.com/files/22629137"

rule chemblDownload:
	input:
		"source/chembl/chembl_30_sqlite.tar.gz"
	output:
		"source/chembl/chembl_30/chembl_30_sqlite/chembl_30.db"
	shell:
		"cd transform/chembl && tar xvzf ../../source/chembl/chembl_30_sqlite.tar.gz -C ../../source/chembl/"

rule chemblDownload_1:
	output:
		"source/chembl/chembl_30_sqlite.tar.gz"
	shell:
		"cd transform/chembl && curl -o ../../source/chembl/chembl_30_sqlite.tar.gz https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/chembl_30_sqlite.tar.gz"

rule chembDrugMechanismExtract:
	input:
		"source/chembl/chembl_30/chembl_30_sqlite/chembl_30.db",
		"src/bmeg/bmeg-dictionary/gdcdictionary/schemas"
	output:
		"output/chembl/chembDrugMechanismExtract.build.protein_drug_association.json.gz"
	shell:
		"sifter run transform/chembl/drug_mechanism.yaml"

rule ensembl:
	input:
		"source/ensembl/Homo_sapiens.GRCh37.87.chr_patch_hapl_scaff.gff3.gz",
		"GRCh37",
		"src/bmeg/bmeg-dictionary/gdcdictionary/schemas",
		"GRCh37",
		"src/bmeg/bmeg-dictionary/gdcdictionary/schemas",
		"GRCh37",
		"src/bmeg/bmeg-dictionary/gdcdictionary/schemas"
	output:
		"output/ensembl/ensembl.exons.exon.json.gz",
		"output/ensembl/ensembl.genes.gene.json.gz",
		"output/ensembl/ensembl.transcripts.transcript.json.gz"
	shell:
		"sifter run transform/ensembl/ensembl_transform.yaml"

rule gene2ensembl:
	input:
		"source/ensembl/gene2ensembl.gz"
	output:
		"source/ensembl/gene2ensembl.translate.link.json.gz"
	shell:
		"sifter run transform/ensembl/gene2ensembl.yaml"

rule hugoDownload:
	output:
		"source/ensembl/Homo_sapiens.GRCh37.85.uniprot.tsv.gz"
	shell:
		"cd transform/ensembl && curl --verbose --progress-bar --ipv4 --connect-timeout 8 --max-time 120 --retry 128 --ftp-ssl --disable-epsv --ftp-pasv ftp://ftp.ensembl.org/pub/grch37/release-96/tsv/homo_sapiens/Homo_sapiens.GRCh37.85.uniprot.tsv.gz --output ../../source/ensembl/Homo_sapiens.GRCh37.85.uniprot.tsv.gz"

rule hugoDownload_1:
	output:
		"source/ensembl/gene2ensembl.gz"
	shell:
		"cd transform/ensembl && curl -o ../../source/ensembl/gene2ensembl.gz https://ftp.ncbi.nih.gov/gene/DATA/gene2ensembl.gz"

rule hugoDownload_2:
	output:
		"source/hugo/hugo.tsv"
	shell:
		"cd transform/ensembl && curl -o ../../source/hugo/hugo.tsv https://www.genenames.org/cgi-bin/download/custom?col=gd_hgnc_id&col=gd_app_sym&col=gd_app_name&col=gd_status&col=gd_prev_sym&col=gd_aliases&col=gd_pub_chrom_map&col=gd_pub_acc_ids&col=gd_pub_refseq_ids&col=md_ensembl_id&col=md_prot_id&status=Approved&status=Entry%20Withdrawn&hgnc_dbtag=on&order_by=gd_app_sym_sort&format=text&submit=submit"

rule hugoMapping:
	input:
		"source/hugo/hugo.tsv"
	output:
		"source/ensembl/uniprotId2ensemblGene.alt.tsv"
	shell:
		"sifter run transform/ensembl/hugo_mapping.yaml"

rule uniprot:
	input:
		"source/ensembl/Homo_sapiens.GRCh37.85.uniprot.tsv.gz",
		"GRCh37"
	output:
		"output/ensembl/uniprot.transform.protein.json.gz"
	shell:
		"sifter run transform/ensembl/uniprot_transform.yaml"

rule download_1:
	output:
		"source/gdc/cases.json"
	shell:
		"cd transform/gdc && ./gdc-scan.py cases ../../source/gdc/cases.json"

rule download_2:
	output:
		"source/gdc/files.json"
	shell:
		"cd transform/gdc && ./gdc-scan.py files ../../source/gdc/files.json"

rule rnaseq:
	input:
		"source/gdc/rna-seq"
	output:
		"output/gdc/rnaseq.rna.debug.json.gz"
	shell:
		"sifter run transform/gdc/rna-expression.yaml"

rule gdc:
	input:
		"source/gdc/cases.json",
		"src/bmeg/bmeg-dictionary/gdcdictionary/schemas",
		"src/bmeg/bmeg-dictionary/gdcdictionary/schemas",
		"src/bmeg/bmeg-dictionary/gdcdictionary/schemas"
	output:
		"output/gdc/gdc.aliquotAlias.table.json.gz",
		"output/gdc/gdc.aliquotObject.aliquot.json.gz",
		"output/gdc/gdc.caseObject.case.json.gz",
		"output/gdc/gdc.sampleObject.sample.json.gz"
	shell:
		"sifter run transform/gdc/transform.yaml"

rule download_3:
	output:
		"source/gdsc/model_list_20191104.csv"
	shell:
		"cd transform/gdsc && curl https://cog.sanger.ac.uk/cmp/download/model_list_20191104.csv -o ../../source/gdsc/model_list_20191104.csv"

rule download_4:
	output:
		"source/go/go.json"
	shell:
		"cd transform/go && curl -o ../../source/go/go.json http://release.geneontology.org/2022-09-19/ontology/go.json"

rule download_5:
	output:
		"source/go/goa_human.gaf.gz"
	shell:
		"cd transform/go && curl -o ../../source/go/goa_human.gaf.gz http://release.geneontology.org/2022-09-19/annotations/goa_human.gaf.gz"

rule go_gaf:
	input:
		"source/go/goa_human.gaf.gz"
	output:
		"output/go/go_gaf.dump.gaf.json.gz"
	shell:
		"sifter run transform/go/go_gaf.yaml"

rule go:
	input:
		"source/go/go.json",
		"src/bmeg/bmeg-dictionary/gdcdictionary/schemas"
	output:
		"output/go/go.edgeOpen.edges.json.gz",
		"output/go/go.nodeOpen.before.json.gz",
		"output/go/go.nodeOpen.term.json.gz"
	shell:
		"sifter run transform/go/go_json.yaml"

rule goReduce:
	input:
		""
	output:
		"output/go/goReduce.edgeReduce.edge_reduce.json.gz"
	shell:
		"sifter run transform/go/go_reduce.yaml"

rule GTEX_Gene_Expression:
	input:
		"source/gtex/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.transpose.gz",
		"src/bmeg/bmeg-dictionary/gdcdictionary/schemas"
	shell:
		"sifter run transform/gtex/gene_transform.yaml"

rule GTEX_Transcript_Expression:
	input:
		"source/gtex/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct.transpose.gz",
		"src/bmeg/bmeg-dictionary/gdcdictionary/schemas"
	shell:
		"sifter run transform/gtex/transcript_transform.yaml"

rule mc3Download:
	output:
		"source/mc3/mc3.v0.2.8.PUBLIC.maf.gz"
	shell:
		"cd transform/mc3 && curl -o ../../source/mc3/mc3.v0.2.8.PUBLIC.maf.gz https://api.gdc.cancer.gov/data/1c8cfe5f-e52d-41ba-94da-f15ea1337efc"

rule mc3:
	input:
		"source/mc3/mc3.v0.2.8.PUBLIC.maf.gz",
		"src/bmeg/bmeg-dictionary/gdcdictionary/schemas",
		"output/gdc/gdc.aliquotAlias.table.json.gz",
		"output/gdc/gdc.aliquotAlias.table.json.gz",
		"output/gdc/gdc.aliquotAlias.table.json.gz"
	output:
		"output/mc3/mc3.callset.callset.json.gz",
		"output/mc3/mc3.transform.raw.json.gz",
		"output/mc3/mc3.variant.allele.json.gz"
	shell:
		"sifter run transform/mc3/transform.yaml"

rule mondo:
	input:
		"source/mondo/mondo.json"
	output:
		"output/mondo/mondo.extractEdges.edges.json.gz",
		"output/mondo/mondo.extractNodes.nodes.json.gz"
	shell:
		"sifter run transform/mondo/mondo.yaml"

rule downloads:
	output:
		"source/msigdb/msigdb_v7.5.1.xml"
	shell:
		"cd transform/msigdb && curl -o ../../source/msigdb/msigdb_v7.5.1.xml https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.5.1/msigdb_v7.5.1.xml"

rule msigdb:
	input:
		"source/msigdb/msigdb_v7.5.1.xml",
		"src/bmeg/bmeg-dictionary/gdcdictionary/schemas",
		"source/ensembl/gene2ensembl.translate.link.json.gz"
	output:
		"output/msigdb/msigdb.transform.raw.json.gz",
		"output/msigdb/msigdb.transform.record.json.gz",
		"output/msigdb/msigdb.object.gene_set.json.gz"
	shell:
		"sifter run transform/msigdb/transform.yaml"

rule :
	input:
		"transform/pathway_commons/target/pc-extract-1.0-jar-with-dependencies.jar",
		"source/pathway_commons/PathwayCommons12.mirtarbase.BIOPAX.owl.gz"
	output:
		"source/pathway_commons/PathwayCommons12.mirtarbase.extSIF",
		"source/pathway_commons/PathwayCommons12.mirtarbase.complex"
	resources:
		mem_mb=15000
	shell:
		"cd transform/pathway_commons && java -Xmx15g -jar target/pc-extract-1.0-jar-with-dependencies.jar ../../source/pathway_commons/PathwayCommons12.mirtarbase.BIOPAX.owl.gz ../../source/pathway_commons/PathwayCommons12.mirtarbase"

rule _1:
	input:
		"transform/pathway_commons/target/pc-extract-1.0-jar-with-dependencies.jar",
		"source/pathway_commons/PathwayCommons12.netpath.BIOPAX.owl.gz"
	output:
		"source/pathway_commons/PathwayCommons12.netpath.extSIF",
		"source/pathway_commons/PathwayCommons12.netpath.complex"
	resources:
		mem_mb=15000
	shell:
		"cd transform/pathway_commons && java -Xmx15g -jar target/pc-extract-1.0-jar-with-dependencies.jar ../../source/pathway_commons/PathwayCommons12.netpath.BIOPAX.owl.gz ../../source/pathway_commons/PathwayCommons12.netpath"

rule _2:
	output:
		"source/pathway_commons/PathwayCommons12.reconx.BIOPAX.owl.gz"
	shell:
		"cd transform/pathway_commons && curl -o ../../source/pathway_commons/PathwayCommons12.reconx.BIOPAX.owl.gz https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.reconx.BIOPAX.owl.gz"

rule _3:
	input:
		"transform/pathway_commons/target/pc-extract-1.0-jar-with-dependencies.jar",
		"source/pathway_commons/PathwayCommons12.bind.BIOPAX.owl.gz"
	output:
		"source/pathway_commons/PathwayCommons12.bind.extSIF",
		"source/pathway_commons/PathwayCommons12.bind.complex"
	resources:
		mem_mb=15000
	shell:
		"cd transform/pathway_commons && java -Xmx15g -jar target/pc-extract-1.0-jar-with-dependencies.jar ../../source/pathway_commons/PathwayCommons12.bind.BIOPAX.owl.gz ../../source/pathway_commons/PathwayCommons12.bind"

rule _4:
	input:
		"transform/pathway_commons/target/pc-extract-1.0-jar-with-dependencies.jar",
		"source/pathway_commons/PathwayCommons12.corum.BIOPAX.owl.gz"
	output:
		"source/pathway_commons/PathwayCommons12.corum.extSIF",
		"source/pathway_commons/PathwayCommons12.corum.complex"
	resources:
		mem_mb=15000
	shell:
		"cd transform/pathway_commons && java -Xmx15g -jar target/pc-extract-1.0-jar-with-dependencies.jar ../../source/pathway_commons/PathwayCommons12.corum.BIOPAX.owl.gz ../../source/pathway_commons/PathwayCommons12.corum"

rule _5:
	input:
		"transform/pathway_commons/target/pc-extract-1.0-jar-with-dependencies.jar",
		"source/pathway_commons/PathwayCommons12.pid.BIOPAX.owl.gz"
	output:
		"source/pathway_commons/PathwayCommons12.pid.extSIF",
		"source/pathway_commons/PathwayCommons12.pid.complex"
	resources:
		mem_mb=15000
	shell:
		"cd transform/pathway_commons && java -Xmx15g -jar target/pc-extract-1.0-jar-with-dependencies.jar ../../source/pathway_commons/PathwayCommons12.pid.BIOPAX.owl.gz ../../source/pathway_commons/PathwayCommons12.pid"

rule _6:
	input:
		"transform/pathway_commons/target/pc-extract-1.0-jar-with-dependencies.jar",
		"source/pathway_commons/PathwayCommons12.psp.BIOPAX.owl.gz"
	output:
		"source/pathway_commons/PathwayCommons12.psp.extSIF",
		"source/pathway_commons/PathwayCommons12.psp.complex"
	resources:
		mem_mb=15000
	shell:
		"cd transform/pathway_commons && java -Xmx15g -jar target/pc-extract-1.0-jar-with-dependencies.jar ../../source/pathway_commons/PathwayCommons12.psp.BIOPAX.owl.gz ../../source/pathway_commons/PathwayCommons12.psp"

rule _7:
	input:
		"transform/pathway_commons/target/pc-extract-1.0-jar-with-dependencies.jar",
		"source/pathway_commons/PathwayCommons12.hprd.BIOPAX.owl.gz"
	output:
		"source/pathway_commons/PathwayCommons12.hprd.extSIF",
		"source/pathway_commons/PathwayCommons12.hprd.complex"
	resources:
		mem_mb=15000
	shell:
		"cd transform/pathway_commons && java -Xmx15g -jar target/pc-extract-1.0-jar-with-dependencies.jar ../../source/pathway_commons/PathwayCommons12.hprd.BIOPAX.owl.gz ../../source/pathway_commons/PathwayCommons12.hprd"

rule _8:
	output:
		"source/pathway_commons/PathwayCommons12.pathbank.BIOPAX.owl.gz"
	shell:
		"cd transform/pathway_commons && curl -o ../../source/pathway_commons/PathwayCommons12.pathbank.BIOPAX.owl.gz https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.pathbank.BIOPAX.owl.gz"

rule _9:
	input:
		"transform/pathway_commons/target/pc-extract-1.0-jar-with-dependencies.jar",
		"source/pathway_commons/PathwayCommons12.msigdb.BIOPAX.owl.gz"
	output:
		"source/pathway_commons/PathwayCommons12.msigdb.extSIF",
		"source/pathway_commons/PathwayCommons12.msigdb.complex"
	resources:
		mem_mb=15000
	shell:
		"cd transform/pathway_commons && java -Xmx15g -jar target/pc-extract-1.0-jar-with-dependencies.jar ../../source/pathway_commons/PathwayCommons12.msigdb.BIOPAX.owl.gz ../../source/pathway_commons/PathwayCommons12.msigdb"

rule _10:
	input:
		"transform/pathway_commons/target/pc-extract-1.0-jar-with-dependencies.jar",
		"source/pathway_commons/PathwayCommons12.reconx.BIOPAX.owl.gz"
	output:
		"source/pathway_commons/PathwayCommons12.reconx.extSIF",
		"source/pathway_commons/PathwayCommons12.reconx.complex"
	resources:
		mem_mb=15000
	shell:
		"cd transform/pathway_commons && java -Xmx15g -jar target/pc-extract-1.0-jar-with-dependencies.jar ../../source/pathway_commons/PathwayCommons12.reconx.BIOPAX.owl.gz ../../source/pathway_commons/PathwayCommons12.reconx"

rule _11:
	output:
		"source/pathway_commons/PathwayCommons12.mirtarbase.BIOPAX.owl.gz"
	shell:
		"cd transform/pathway_commons && curl -o ../../source/pathway_commons/PathwayCommons12.mirtarbase.BIOPAX.owl.gz https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.mirtarbase.BIOPAX.owl.gz"

rule _12:
	output:
		"source/pathway_commons/PathwayCommons12.pid.BIOPAX.owl.gz"
	shell:
		"cd transform/pathway_commons && curl -o ../../source/pathway_commons/PathwayCommons12.pid.BIOPAX.owl.gz https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.pid.BIOPAX.owl.gz"

rule _13:
	input:
		"transform/pathway_commons/target/pc-extract-1.0-jar-with-dependencies.jar",
		"source/pathway_commons/PathwayCommons12.inoh.BIOPAX.owl.gz"
	output:
		"source/pathway_commons/PathwayCommons12.inoh.extSIF",
		"source/pathway_commons/PathwayCommons12.inoh.complex"
	resources:
		mem_mb=15000
	shell:
		"cd transform/pathway_commons && java -Xmx15g -jar target/pc-extract-1.0-jar-with-dependencies.jar ../../source/pathway_commons/PathwayCommons12.inoh.BIOPAX.owl.gz ../../source/pathway_commons/PathwayCommons12.inoh"

rule _14:
	output:
		"source/pathway_commons/PathwayCommons12.humancyc.BIOPAX.owl.gz"
	shell:
		"cd transform/pathway_commons && curl -o ../../source/pathway_commons/PathwayCommons12.humancyc.BIOPAX.owl.gz https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.humancyc.BIOPAX.owl.gz"

rule _15:
	input:
		"transform/pathway_commons/target/pc-extract-1.0-jar-with-dependencies.jar",
		"source/pathway_commons/PathwayCommons12.ctd.BIOPAX.owl.gz"
	output:
		"source/pathway_commons/PathwayCommons12.ctd.extSIF",
		"source/pathway_commons/PathwayCommons12.ctd.complex"
	resources:
		mem_mb=15000
	shell:
		"cd transform/pathway_commons && java -Xmx15g -jar target/pc-extract-1.0-jar-with-dependencies.jar ../../source/pathway_commons/PathwayCommons12.ctd.BIOPAX.owl.gz ../../source/pathway_commons/PathwayCommons12.ctd"

rule _16:
	input:
		"transform/pathway_commons/target/pc-extract-1.0-jar-with-dependencies.jar",
		"source/pathway_commons/PathwayCommons12.humancyc.BIOPAX.owl.gz"
	output:
		"source/pathway_commons/PathwayCommons12.humancyc.extSIF",
		"source/pathway_commons/PathwayCommons12.humancyc.complex"
	resources:
		mem_mb=15000
	shell:
		"cd transform/pathway_commons && java -Xmx15g -jar target/pc-extract-1.0-jar-with-dependencies.jar ../../source/pathway_commons/PathwayCommons12.humancyc.BIOPAX.owl.gz ../../source/pathway_commons/PathwayCommons12.humancyc"

rule _17:
	output:
		"source/pathway_commons/PathwayCommons12.bind.BIOPAX.owl.gz"
	shell:
		"cd transform/pathway_commons && curl -o ../../source/pathway_commons/PathwayCommons12.bind.BIOPAX.owl.gz https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.bind.BIOPAX.owl.gz"

rule _18:
	output:
		"source/pathway_commons/PathwayCommons12.biogrid.BIOPAX.owl.gz"
	shell:
		"cd transform/pathway_commons && curl -o ../../source/pathway_commons/PathwayCommons12.biogrid.BIOPAX.owl.gz https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.biogrid.BIOPAX.owl.gz"

rule _19:
	input:
		"transform/pathway_commons/target/pc-extract-1.0-jar-with-dependencies.jar",
		"source/pathway_commons/PathwayCommons12.reactome.BIOPAX.owl.gz"
	output:
		"source/pathway_commons/PathwayCommons12.reactome.extSIF",
		"source/pathway_commons/PathwayCommons12.reactome.complex"
	resources:
		mem_mb=15000
	shell:
		"cd transform/pathway_commons && java -Xmx15g -jar target/pc-extract-1.0-jar-with-dependencies.jar ../../source/pathway_commons/PathwayCommons12.reactome.BIOPAX.owl.gz ../../source/pathway_commons/PathwayCommons12.reactome"

rule _20:
	output:
		"source/pathway_commons/PathwayCommons12.ctd.BIOPAX.owl.gz"
	shell:
		"cd transform/pathway_commons && curl -o ../../source/pathway_commons/PathwayCommons12.ctd.BIOPAX.owl.gz https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.ctd.BIOPAX.owl.gz"

rule _21:
	input:
		"transform/pathway_commons/target/pc-extract-1.0-jar-with-dependencies.jar",
		"source/pathway_commons/PathwayCommons12.drugbank.BIOPAX.owl.gz"
	output:
		"source/pathway_commons/PathwayCommons12.drugbank.extSIF",
		"source/pathway_commons/PathwayCommons12.drugbank.complex"
	resources:
		mem_mb=15000
	shell:
		"cd transform/pathway_commons && java -Xmx15g -jar target/pc-extract-1.0-jar-with-dependencies.jar ../../source/pathway_commons/PathwayCommons12.drugbank.BIOPAX.owl.gz ../../source/pathway_commons/PathwayCommons12.drugbank"

rule _22:
	output:
		"source/pathway_commons/PathwayCommons12.kegg.BIOPAX.owl.gz"
	shell:
		"cd transform/pathway_commons && curl -o ../../source/pathway_commons/PathwayCommons12.kegg.BIOPAX.owl.gz https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.kegg.BIOPAX.owl.gz"

rule _23:
	output:
		"source/pathway_commons/PathwayCommons12.dip.BIOPAX.owl.gz"
	shell:
		"cd transform/pathway_commons && curl -o ../../source/pathway_commons/PathwayCommons12.dip.BIOPAX.owl.gz https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.dip.BIOPAX.owl.gz"

rule _24:
	output:
		"source/pathway_commons/PathwayCommons12.inoh.BIOPAX.owl.gz"
	shell:
		"cd transform/pathway_commons && curl -o ../../source/pathway_commons/PathwayCommons12.inoh.BIOPAX.owl.gz https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.inoh.BIOPAX.owl.gz"

rule _25:
	input:
		"transform/pathway_commons/target/pc-extract-1.0-jar-with-dependencies.jar",
		"source/pathway_commons/PathwayCommons12.biogrid.BIOPAX.owl.gz"
	output:
		"source/pathway_commons/PathwayCommons12.biogrid.extSIF",
		"source/pathway_commons/PathwayCommons12.biogrid.complex"
	resources:
		mem_mb=15000
	shell:
		"cd transform/pathway_commons && java -Xmx15g -jar target/pc-extract-1.0-jar-with-dependencies.jar ../../source/pathway_commons/PathwayCommons12.biogrid.BIOPAX.owl.gz ../../source/pathway_commons/PathwayCommons12.biogrid"

rule _26:
	input:
		"transform/pathway_commons/target/pc-extract-1.0-jar-with-dependencies.jar",
		"source/pathway_commons/PathwayCommons12.pathbank.BIOPAX.owl.gz"
	output:
		"source/pathway_commons/PathwayCommons12.pathbank.extSIF",
		"source/pathway_commons/PathwayCommons12.pathbank.complex"
	resources:
		mem_mb=15000
	shell:
		"cd transform/pathway_commons && java -Xmx15g -jar target/pc-extract-1.0-jar-with-dependencies.jar ../../source/pathway_commons/PathwayCommons12.pathbank.BIOPAX.owl.gz ../../source/pathway_commons/PathwayCommons12.pathbank"

rule _27:
	input:
		"transform/pathway_commons/target/pc-extract-1.0-jar-with-dependencies.jar",
		"source/pathway_commons/PathwayCommons12.innatedb.BIOPAX.owl.gz"
	output:
		"source/pathway_commons/PathwayCommons12.innatedb.extSIF",
		"source/pathway_commons/PathwayCommons12.innatedb.complex"
	resources:
		mem_mb=15000
	shell:
		"cd transform/pathway_commons && java -Xmx15g -jar target/pc-extract-1.0-jar-with-dependencies.jar ../../source/pathway_commons/PathwayCommons12.innatedb.BIOPAX.owl.gz ../../source/pathway_commons/PathwayCommons12.innatedb"

rule _28:
	output:
		"source/pathway_commons/PathwayCommons12.msigdb.BIOPAX.owl.gz"
	shell:
		"cd transform/pathway_commons && curl -o ../../source/pathway_commons/PathwayCommons12.msigdb.BIOPAX.owl.gz https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.msigdb.BIOPAX.owl.gz"

rule _29:
	input:
		"transform/pathway_commons/target/pc-extract-1.0-jar-with-dependencies.jar",
		"source/pathway_commons/PathwayCommons12.kegg.BIOPAX.owl.gz"
	output:
		"source/pathway_commons/PathwayCommons12.kegg.extSIF",
		"source/pathway_commons/PathwayCommons12.kegg.complex"
	resources:
		mem_mb=15000
	shell:
		"cd transform/pathway_commons && java -Xmx15g -jar target/pc-extract-1.0-jar-with-dependencies.jar ../../source/pathway_commons/PathwayCommons12.kegg.BIOPAX.owl.gz ../../source/pathway_commons/PathwayCommons12.kegg"

rule _30:
	output:
		"source/pathway_commons/PathwayCommons12.psp.BIOPAX.owl.gz"
	shell:
		"cd transform/pathway_commons && curl -o ../../source/pathway_commons/PathwayCommons12.psp.BIOPAX.owl.gz https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.psp.BIOPAX.owl.gz"

rule _31:
	input:
		"transform/pathway_commons/target/pc-extract-1.0-jar-with-dependencies.jar",
		"source/pathway_commons/PathwayCommons12.panther.BIOPAX.owl.gz"
	output:
		"source/pathway_commons/PathwayCommons12.panther.extSIF",
		"source/pathway_commons/PathwayCommons12.panther.complex"
	resources:
		mem_mb=15000
	shell:
		"cd transform/pathway_commons && java -Xmx15g -jar target/pc-extract-1.0-jar-with-dependencies.jar ../../source/pathway_commons/PathwayCommons12.panther.BIOPAX.owl.gz ../../source/pathway_commons/PathwayCommons12.panther"

rule _32:
	output:
		"source/pathway_commons/PathwayCommons12.hprd.BIOPAX.owl.gz"
	shell:
		"cd transform/pathway_commons && curl -o ../../source/pathway_commons/PathwayCommons12.hprd.BIOPAX.owl.gz https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.hprd.BIOPAX.owl.gz"

rule _33:
	output:
		"source/pathway_commons/PathwayCommons12.reactome.BIOPAX.owl.gz"
	shell:
		"cd transform/pathway_commons && curl -o ../../source/pathway_commons/PathwayCommons12.reactome.BIOPAX.owl.gz https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.reactome.BIOPAX.owl.gz"

rule _34:
	input:
		"transform/pathway_commons/target/pc-extract-1.0-jar-with-dependencies.jar",
		"source/pathway_commons/PathwayCommons12.dip.BIOPAX.owl.gz"
	output:
		"source/pathway_commons/PathwayCommons12.dip.extSIF",
		"source/pathway_commons/PathwayCommons12.dip.complex"
	resources:
		mem_mb=15000
	shell:
		"cd transform/pathway_commons && java -Xmx15g -jar target/pc-extract-1.0-jar-with-dependencies.jar ../../source/pathway_commons/PathwayCommons12.dip.BIOPAX.owl.gz ../../source/pathway_commons/PathwayCommons12.dip"

rule _35:
	output:
		"source/pathway_commons/PathwayCommons12.drugbank.BIOPAX.owl.gz"
	shell:
		"cd transform/pathway_commons && curl -o ../../source/pathway_commons/PathwayCommons12.drugbank.BIOPAX.owl.gz https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.drugbank.BIOPAX.owl.gz"

rule _36:
	output:
		"source/pathway_commons/PathwayCommons12.corum.BIOPAX.owl.gz"
	shell:
		"cd transform/pathway_commons && curl -o ../../source/pathway_commons/PathwayCommons12.corum.BIOPAX.owl.gz https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.corum.BIOPAX.owl.gz"

rule _37:
	output:
		"transform/pathway_commons/target/pc-extract-1.0-jar-with-dependencies.jar"
	shell:
		"cd transform/pathway_commons && docker run  -u `id -u` -v `pwd`:`pwd` -v `pwd`/maven-repo:/root/.m2 -w `pwd` -ti maven:3-openjdk-11 mvn package"

rule _38:
	output:
		"source/pathway_commons/PathwayCommons12.netpath.BIOPAX.owl.gz"
	shell:
		"cd transform/pathway_commons && curl -o ../../source/pathway_commons/PathwayCommons12.netpath.BIOPAX.owl.gz https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.netpath.BIOPAX.owl.gz"

rule _39:
	output:
		"source/pathway_commons/PathwayCommons12.innatedb.BIOPAX.owl.gz"
	shell:
		"cd transform/pathway_commons && curl -o ../../source/pathway_commons/PathwayCommons12.innatedb.BIOPAX.owl.gz https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.innatedb.BIOPAX.owl.gz"

rule _40:
	output:
		"source/pathway_commons/PathwayCommons12.panther.BIOPAX.owl.gz"
	shell:
		"cd transform/pathway_commons && curl -o ../../source/pathway_commons/PathwayCommons12.panther.BIOPAX.owl.gz https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.panther.BIOPAX.owl.gz"

rule pfam_family:
	input:
		"source/pfam/xmls",
		"src/bmeg/bmeg-dictionary/gdcdictionary/schemas"
	output:
		"output/pfam/pfam_family.transform.pfam.json.gz",
		"output/pfam/pfam_family.transform.debug.json.gz"
	shell:
		"sifter run transform/pfam/transform.yaml"

rule pharmacodb:
	input:
		"source/pharmacodb/pharmacodb-1.1.1.sql"
	output:
		"output/pharmacodb/tissues.tsv.gz",
		"output/pharmacodb/experiments.tsv.gz",
		"output/pharmacodb/profiles.tsv.gz",
		"output/pharmacodb/drugs.tsv.gz",
		"output/pharmacodb/source_tissue_names.tsv.gz",
		"output/pharmacodb/cells.tsv.gz",
		"output/pharmacodb/dose_responses.tsv.gz",
		"output/pharmacodb/datasets.tsv.gz",
		"output/pharmacodb/source_statistics.tsv.gz",
		"output/pharmacodb/source_cell_names.tsv.gz",
		"output/pharmacodb/source_drug_names.tsv.gz",
		"output/pharmacodb/cellosaurus.tsv.gz",
		"output/pharmacodb/cell_tissues.tsv.gz",
		"output/pharmacodb/sources.tsv.gz",
		"output/pharmacodb/dataset_cells.tsv.gz",
		"output/pharmacodb/drug_annots.tsv.gz"
	shell:
		"sifter run transform/pharmacodb/database.yaml"

rule pharmacodbDownload:
	output:
		"source/pharmacodb/pharmacodb-1.1.1.sql"
	shell:
		"cd transform/pharmacodb && curl -o ../../source/pharmacodb/pharmacodb-1.1.1.sql https://zenodo.org/record/1143645/files/pharmacodb-1.1.1.sql"

rule pharmacodb_profiles:
	input:
		"output/pharmacodb/profiles.tsv.gz",
		"src/bmeg/bmeg-dictionary/gdcdictionary/schemas",
		"src/bmeg/bmeg-dictionary/gdcdictionary/schemas",
		"src/bmeg/bmeg-dictionary/gdcdictionary/schemas",
		"src/bmeg/bmeg-dictionary/gdcdictionary/schemas",
		"src/bmeg/bmeg-dictionary/gdcdictionary/schemas",
		"output/pharmacodb/dr_reduce.curveReduce.dose_response_curve.json.gz",
		"output/pharmacodb/experiments.tsv.gz",
		"output/pharmacodb/drugs.tsv.gz",
		"output/pharmacodb/datasets.tsv.gz",
		"output/pharmacodb/cells.tsv.gz",
		"source/ccle/cellline_id_lookup.tsv"
	output:
		"output/pharmacodb/pharmacodb_profiles.cellCase.case.json.gz",
		"output/pharmacodb/pharmacodb_profiles.cellDistinct.checkpoint.json.gz",
		"output/pharmacodb/pharmacodb_profiles.cellProject.project.json.gz",
		"output/pharmacodb/pharmacodb_profiles.cellSample.sample.json.gz",
		"output/pharmacodb/pharmacodb_profiles.drObject.drug_response.json.gz",
		"output/pharmacodb/pharmacodb_profiles.cellAliquot.aliquot.json.gz"
	shell:
		"sifter run transform/pharmacodb/profile.yaml"

rule dr_reduce:
	input:
		"output/pharmacodb/dose_responses.tsv.gz"
	output:
		"output/pharmacodb/dr_reduce.curveReduce.dose_response_curve.json.gz"
	shell:
		"sifter run transform/pharmacodb/reduce_dr_data.yaml"

rule pubmed:
	input:
		"source/pubmed/baseline",
		"src/bmeg/bmeg-dictionary/gdcdictionary/schemas"
	output:
		"output/pubmed/pubmed.transform.publication.json.gz"
	shell:
		"sifter run transform/pubmed/transform.yaml"

rule chembDrugMechanismExtract_1:
	input:
		"output/chembl/chembDrugMechanismExtract.build.protein_drug_association.json.gz",
		"src/bmeg/bmeg-dictionary/gdcdictionary/schemas"
	output:
		"graph/chembDrugMechanismExtract.build-graph.edge.json.gz",
		"graph/chembDrugMechanismExtract.build-graph.vertex.json.gz"
	shell:
		"sifter run graph-build/chembDrugMechanismExtract.yaml"

rule ensembl_1:
	input:
		"output/ensembl/ensembl.transcripts.transcript.json.gz",
		"output/ensembl/ensembl.exons.exon.json.gz",
		"output/ensembl/ensembl.genes.gene.json.gz",
		"src/bmeg/bmeg-dictionary/gdcdictionary/schemas",
		"src/bmeg/bmeg-dictionary/gdcdictionary/schemas",
		"src/bmeg/bmeg-dictionary/gdcdictionary/schemas"
	output:
		"graph/ensembl.genes-graph.edge.json.gz",
		"graph/ensembl.transcripts-graph.vertex.json.gz",
		"graph/ensembl.transcripts-graph.edge.json.gz",
		"graph/ensembl.exons-graph.vertex.json.gz",
		"graph/ensembl.exons-graph.edge.json.gz",
		"graph/ensembl.genes-graph.vertex.json.gz"
	shell:
		"sifter run graph-build/ensembl.yaml"

rule gdc_1:
	input:
		"output/gdc/gdc.aliquotObject.aliquot.json.gz",
		"output/gdc/gdc.caseObject.case.json.gz",
		"output/gdc/gdc.sampleObject.sample.json.gz",
		"src/bmeg/bmeg-dictionary/gdcdictionary/schemas",
		"src/bmeg/bmeg-dictionary/gdcdictionary/schemas",
		"src/bmeg/bmeg-dictionary/gdcdictionary/schemas"
	output:
		"graph/gdc.aliquotObject-graph.vertex.json.gz",
		"graph/gdc.aliquotObject-graph.edge.json.gz",
		"graph/gdc.caseObject-graph.vertex.json.gz",
		"graph/gdc.caseObject-graph.edge.json.gz",
		"graph/gdc.sampleObject-graph.vertex.json.gz",
		"graph/gdc.sampleObject-graph.edge.json.gz"
	shell:
		"sifter run graph-build/gdc.yaml"

rule go_1:
	input:
		"output/go/go.nodeOpen.term.json.gz",
		"src/bmeg/bmeg-dictionary/gdcdictionary/schemas"
	output:
		"graph/go.nodeOpen-graph.edge.json.gz",
		"graph/go.nodeOpen-graph.vertex.json.gz"
	shell:
		"sifter run graph-build/go.yaml"

rule mc3_1:
	input:
		"output/mc3/mc3.callset.callset.json.gz",
		"src/bmeg/bmeg-dictionary/gdcdictionary/schemas"
	output:
		"graph/mc3.callset-graph.edge.json.gz",
		"graph/mc3.callset-graph.vertex.json.gz"
	shell:
		"sifter run graph-build/mc3.yaml"

rule msigdb_1:
	input:
		"output/msigdb/msigdb.object.gene_set.json.gz",
		"src/bmeg/bmeg-dictionary/gdcdictionary/schemas"
	output:
		"graph/msigdb.object-graph.vertex.json.gz",
		"graph/msigdb.object-graph.edge.json.gz"
	shell:
		"sifter run graph-build/msigdb.yaml"

rule pfam_family_1:
	input:
		"output/pfam/pfam_family.transform.pfam.json.gz",
		"src/bmeg/bmeg-dictionary/gdcdictionary/schemas"
	output:
		"graph/pfam_family.transform-graph.vertex.json.gz",
		"graph/pfam_family.transform-graph.edge.json.gz"
	shell:
		"sifter run graph-build/pfam_family.yaml"

rule pharmacodb_profiles_1:
	input:
		"output/pharmacodb/pharmacodb_profiles.cellAliquot.aliquot.json.gz",
		"output/pharmacodb/pharmacodb_profiles.cellCase.case.json.gz",
		"output/pharmacodb/pharmacodb_profiles.cellProject.project.json.gz",
		"output/pharmacodb/pharmacodb_profiles.cellSample.sample.json.gz",
		"output/pharmacodb/pharmacodb_profiles.drObject.drug_response.json.gz",
		"src/bmeg/bmeg-dictionary/gdcdictionary/schemas",
		"src/bmeg/bmeg-dictionary/gdcdictionary/schemas",
		"src/bmeg/bmeg-dictionary/gdcdictionary/schemas",
		"src/bmeg/bmeg-dictionary/gdcdictionary/schemas",
		"src/bmeg/bmeg-dictionary/gdcdictionary/schemas"
	output:
		"graph/pharmacodb_profiles.drObject-graph.vertex.json.gz",
		"graph/pharmacodb_profiles.drObject-graph.edge.json.gz",
		"graph/pharmacodb_profiles.cellCase-graph.vertex.json.gz",
		"graph/pharmacodb_profiles.cellProject-graph.edge.json.gz",
		"graph/pharmacodb_profiles.cellSample-graph.vertex.json.gz",
		"graph/pharmacodb_profiles.cellAliquot-graph.vertex.json.gz",
		"graph/pharmacodb_profiles.cellAliquot-graph.edge.json.gz",
		"graph/pharmacodb_profiles.cellCase-graph.edge.json.gz",
		"graph/pharmacodb_profiles.cellProject-graph.vertex.json.gz",
		"graph/pharmacodb_profiles.cellSample-graph.edge.json.gz"
	shell:
		"sifter run graph-build/pharmacodb_profiles.yaml"

rule pubmed_1:
	input:
		"output/pubmed/pubmed.transform.publication.json.gz",
		"src/bmeg/bmeg-dictionary/gdcdictionary/schemas"
	output:
		"graph/pubmed.transform-graph.vertex.json.gz",
		"graph/pubmed.transform-graph.edge.json.gz"
	shell:
		"sifter run graph-build/pubmed.yaml"



