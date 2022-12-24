


rule all:
	input:
		"output/gdc/open-maf/gdc-mafs.scan.variant.json.gz",
		"source/pathway_commons/PathwayCommons12.hprd.complex",
		"output/depmap/depmap-expression.values.expression.json.gz",
		"output/gdc/gdc.sampleObject.sample.json.gz",
		"output/pharmacodb/pharmacodb_profiles.cellProject.project.json.gz",
		"output/depmap/depmap-mafs.callsets.callset.json.gz",
		"source/pathway_commons/PathwayCommons12.humancyc.complex",
		"source/pharmacodb/tables/source_cell_names.tsv.gz",
		"source/pathway_commons/PathwayCommons12.corum.extSIF",
		"source/pathway_commons/PathwayCommons12.innatedb.complex",
		"output/depmap/depmap-cases.cases.case.json.gz",
		"output/go/goJson.edgeOpen.edges.json.gz",
		"source/gtex/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz",
		"source/pathway_commons/PathwayCommons12.hprd.extSIF",
		"source/ccle/sample_info.csv",
		"source/pathway_commons/PathwayCommons12.psp.complex",
		"source/pathway_commons/PathwayCommons12.inoh.complex",
		"output/pharmacodb/pharmacodb_profiles.cellSample.sample.json.gz",
		"output/pathway_commons/pathway_commons.interactionMap.interaction.json.gz",
		"output/ensembl/uniprot.transform.protein.json.gz",
		"source/gtex/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt",
		"source/pathway_commons/PathwayCommons12.panther.complex",
		"source/pathway_commons/PathwayCommons12.reactome.complex",
		"output/ensembl/ensembl.transcripts.transcript.json.gz",
		"output/msigdb/msigdb.transform.raw.json.gz",
		"source/pathway_commons/PathwayCommons12.biogrid.extSIF",
		"source/pathway_commons/PathwayCommons12.ctd.complex",
		"source/pharmacodb/tables/tissues.tsv.gz",
		"output/pharmacodb/pharmacodb_profiles.drObject.drug_response.json.gz",
		"source/pathway_commons/PathwayCommons12.psp.extSIF",
		"source/pathway_commons/PathwayCommons12.pid.extSIF",
		"source/pharmacodb/tables/source_statistics.tsv.gz",
		"source/pathway_commons/PathwayCommons12.biogrid.complex",
		"source/pathway_commons/PathwayCommons12.netpath.complex",
		"source/gdsc/model_list_20191104.csv",
		"source/pathway_commons/PathwayCommons12.panther.extSIF",
		"output/gdc/gdc.aliquotAlias.table.json.gz",
		"output/mondo/mondo.extractNodes.nodes.json.gz",
		"output/msigdb/msigdb.object.gene_set.json.gz",
		"source/pathway_commons/PathwayCommons12.bind.extSIF",
		"source/pathway_commons/PathwayCommons12.msigdb.extSIF",
		"source/pathway_commons/PathwayCommons12.kegg.complex",
		"source/pharmacodb/tables/source_tissue_names.tsv.gz",
		"source/pathway_commons/PathwayCommons12.bind.complex",
		"source/pathway_commons/PathwayCommons12.kegg.extSIF",
		"source/pathway_commons/PathwayCommons12.pathbank.extSIF",
		"source/pathway_commons/PathwayCommons12.netpath.extSIF",
		"source/pathway_commons/PathwayCommons12.drugbank.extSIF",
		"source/pathway_commons/PathwayCommons12.humancyc.extSIF",
		"output/pathway_commons/pathway_commons.complexBundle.complex.json.gz",
		"output/ensembl/ensembl.genes.gene.json.gz",
		"source/pathway_commons/PathwayCommons12.drugbank.complex",
		"source/ccle/cellline_id_lookup.tsv",
		"output/go/goJson.nodeOpen.term.json.gz",
		"source/pathway_commons/PathwayCommons12.inoh.extSIF",
		"source/pathway_commons/PathwayCommons12.ctd.extSIF",
		"output/go/go_gaf.dump.gaf.json.gz",
		"source/pathway_commons/PathwayCommons12.corum.complex",
		"output/pharmacodb/pharmacodb_profiles.cellDistinct.checkpoint.json.gz",
		"output/gdc/gdc.aliquotObject.aliquot.json.gz",
		"output/gdc/gdc.caseObject.case.json.gz",
		"source/pathway_commons/PathwayCommons12.innatedb.extSIF",
		"source/pharmacodb/tables/cellosaurus.tsv.gz",
		"source/pharmacodb/tables/source_drug_names.tsv.gz",
		"output/pubmed/pubmed.transform.publication.json.gz",
		"source/pathway_commons/PathwayCommons12.reactome.extSIF",
		"output/ensembl/ensembl.exons.exon.json.gz",
		"source/pathway_commons/PathwayCommons12.dip.complex",
		"source/pharmacodb/tables/sources.tsv.gz",
		"source/pharmacodb/tables/drug_annots.tsv.gz",
		"output/chembl/chembDrugMechanismExtract.build.protein_drug_association.json.gz",
		"source/gdc/files.json",
		"source/gtex/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt",
		"output/msigdb/msigdb.transform.record.json.gz",
		"source/pathway_commons/PathwayCommons12.msigdb.complex",
		"source/pathway_commons/PathwayCommons12.mirtarbase.complex",
		"output/pharmacodb/pharmacodb_profiles.cellCase.case.json.gz",
		"output/depmap/depmap-mafs.variants.variants.json.gz",
		"output/gdc/open-maf/gdc-mafs.somaticCallsets.callset.json.gz",
		"source/pathway_commons/PathwayCommons12.reconx.extSIF",
		"source/pharmacodb/tables/dataset_cells.tsv.gz",
		"source/pharmacodb/tables/cell_tissues.tsv.gz",
		"output/pharmacodb/pharmacodb_profiles.cellAliquot.aliquot.json.gz",
		"source/ensembl/uniprotId2ensemblGene.alt.tsv",
		"output/gdc/rnaseq.rna.gene_rnaseq.json.gz",
		"output/go/goJson.nodeOpen.before.json.gz",
		"source/pathway_commons/PathwayCommons12.reconx.complex",
		"source/pathway_commons/PathwayCommons12.dip.extSIF",
		"source/pathway_commons/PathwayCommons12.mirtarbase.extSIF",
		"output/mondo/mondo.extractEdges.edges.json.gz",
		"source/pathway_commons/PathwayCommons12.pid.complex",
		"output/go/goReduce.edgeReduce.edge_reduce.json.gz",
		"source/pathway_commons/PathwayCommons12.pathbank.complex"

rule cell_line_names:
	input:
		"output/pharmacodb/cellosaurus.tsv.gz",
		"transform/ccle/cell_lines.yaml"
	output:
		"source/ccle/cellline_id_lookup.tsv"
	shell:
		"sifter run transform/ccle/cell_lines.yaml"

rule download_sample_info:
	output:
		"source/ccle/sample_info.csv"
	shell:
		"cd transform/ccle && curl -L -o ../../source/ccle/sample_info.csv https://ndownloader.figshare.com/files/22629137"

rule chemblDownload_curl:
	output:
		"source/chembl/chembl_30_sqlite.tar.gz"
	shell:
		"cd transform/chembl && curl -o ../../source/chembl/chembl_30_sqlite.tar.gz https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/chembl_30_sqlite.tar.gz"

rule chemblDownload_tar:
	input:
		"source/chembl/chembl_30_sqlite.tar.gz"
	output:
		"source/chembl/chembl_30/chembl_30_sqlite/chembl_30.db"
	shell:
		"cd transform/chembl && tar xvzf ../../source/chembl/chembl_30_sqlite.tar.gz -C ../../source/chembl/"

rule chembDrugMechanismExtract:
	input:
		"source/chembl/chembl_30/chembl_30_sqlite/chembl_30.db",
		"schema",
		"transform/chembl/drug_mechanism.yaml"
	output:
		"output/chembl/chembDrugMechanismExtract.build.protein_drug_association.json.gz"
	shell:
		"sifter run transform/chembl/drug_mechanism.yaml"

rule depmap_cases:
	input:
		"source/depmap/Model.csv",
		"schema",
		"transform/depmap/cases.yaml"
	output:
		"output/depmap/depmap-cases.cases.case.json.gz"
	shell:
		"sifter run transform/depmap/cases.yaml"

rule depmap_expression:
	input:
		"source/depmap/OmicsExpressionProteinCodingGenesTPMLogp1.csv",
		"source/ensembl/gene2ensembl.translate.link.json.gz",
		"schema",
		"transform/depmap/expression.yaml"
	output:
		"output/depmap/depmap-expression.values.expression.json.gz"
	shell:
		"sifter run transform/depmap/expression.yaml"

rule depmap_mafs:
	input:
		"source/depmap/OmicsSomaticMutations.csv",
		"schema",
		"schema",
		"transform/depmap/mutations.yaml"
	output:
		"output/depmap/depmap-mafs.callsets.callset.json.gz",
		"output/depmap/depmap-mafs.variants.variants.json.gz"
	shell:
		"sifter run transform/depmap/mutations.yaml"

rule ensembl:
	input:
		"source/ensembl/Homo_sapiens.GRCh37.87.chr_patch_hapl_scaff.gff3.gz",
		"schema",
		"schema",
		"schema",
		"transform/ensembl/ensembl_transform.yaml"
	output:
		"output/ensembl/ensembl.transcripts.transcript.json.gz",
		"output/ensembl/ensembl.exons.exon.json.gz",
		"output/ensembl/ensembl.genes.gene.json.gz"
	shell:
		"sifter run transform/ensembl/ensembl_transform.yaml"

rule gene2ensembl:
	input:
		"source/ensembl/gene2ensembl.gz",
		"transform/ensembl/gene2ensembl.yaml"
	output:
		"source/ensembl/gene2ensembl.translate.link.json.gz"
	shell:
		"sifter run transform/ensembl/gene2ensembl.yaml"

rule hugoDownload_gene2entrez:
	output:
		"source/ensembl/gene2ensembl.gz"
	shell:
		"cd transform/ensembl && curl -o ../../source/ensembl/gene2ensembl.gz https://ftp.ncbi.nih.gov/gene/DATA/gene2ensembl.gz"

rule hugoDownload_hugoDownload:
	output:
		"source/hugo/hugo.tsv"
	shell:
		"cd transform/ensembl && curl -o ../../source/hugo/hugo.tsv https://www.genenames.org/cgi-bin/download/custom?col=gd_hgnc_id&col=gd_app_sym&col=gd_app_name&col=gd_status&col=gd_prev_sym&col=gd_aliases&col=gd_pub_chrom_map&col=gd_pub_acc_ids&col=gd_pub_refseq_ids&col=md_ensembl_id&col=md_prot_id&status=Approved&status=Entry%20Withdrawn&hgnc_dbtag=on&order_by=gd_app_sym_sort&format=text&submit=submit"

rule hugoDownload_uniprotDownload:
	output:
		"source/ensembl/Homo_sapiens.GRCh37.85.uniprot.tsv.gz"
	shell:
		"cd transform/ensembl && curl --verbose --progress-bar --ipv4 --connect-timeout 8 --max-time 120 --retry 128 --ftp-ssl --disable-epsv --ftp-pasv ftp://ftp.ensembl.org/pub/grch37/release-96/tsv/homo_sapiens/Homo_sapiens.GRCh37.85.uniprot.tsv.gz --output ../../source/ensembl/Homo_sapiens.GRCh37.85.uniprot.tsv.gz"

rule hugoMapping:
	input:
		"source/hugo/hugo.tsv",
		"transform/ensembl/hugo_mapping.yaml"
	output:
		"source/ensembl/uniprotId2ensemblGene.alt.tsv"
	shell:
		"sifter run transform/ensembl/hugo_mapping.yaml"

rule uniprot:
	input:
		"source/ensembl/Homo_sapiens.GRCh37.85.uniprot.tsv.gz",
		"schemas",
		"transform/ensembl/uniprot_transform.yaml"
	output:
		"output/ensembl/uniprot.transform.protein.json.gz"
	shell:
		"sifter run transform/ensembl/uniprot_transform.yaml"

rule download_cases:
	output:
		"source/gdc/cases.json"
	shell:
		"cd transform/gdc && ./gdc-scan.py cases ../../source/gdc/cases.json"

rule download_files:
	output:
		"source/gdc/files.json"
	shell:
		"cd transform/gdc && ./gdc-scan.py files ../../source/gdc/files.json"

rule gdc_mafs:
	input:
		"source/gdc/open-maf",
		"schema",
		"schema",
		"transform/gdc/maf-files.yaml"
	output:
		"output/gdc/open-maf/gdc-mafs.scan.variant.json.gz",
		"output/gdc/open-maf/gdc-mafs.somaticCallsets.callset.json.gz"
	shell:
		"sifter run transform/gdc/maf-files.yaml"

rule rnaseq:
	input:
		"source/gdc/rna-seq",
		"schema",
		"transform/gdc/rna-expression.yaml"
	output:
		"output/gdc/rnaseq.rna.gene_rnaseq.json.gz"
	shell:
		"sifter run transform/gdc/rna-expression.yaml"

rule gdc:
	input:
		"source/gdc/cases.json",
		"schema",
		"schema",
		"schema",
		"transform/gdc/transform.yaml"
	output:
		"output/gdc/gdc.aliquotObject.aliquot.json.gz",
		"output/gdc/gdc.caseObject.case.json.gz",
		"output/gdc/gdc.sampleObject.sample.json.gz",
		"output/gdc/gdc.aliquotAlias.table.json.gz"
	shell:
		"sifter run transform/gdc/transform.yaml"

rule download_sampleInfo:
	output:
		"source/gdsc/model_list_20191104.csv"
	shell:
		"cd transform/gdsc && curl https://cog.sanger.ac.uk/cmp/download/model_list_20191104.csv -o ../../source/gdsc/model_list_20191104.csv"

rule download_downloadGAF:
	output:
		"source/go/goa_human.gaf.gz"
	shell:
		"cd transform/go && curl -o ../../source/go/goa_human.gaf.gz http://release.geneontology.org/2022-09-19/annotations/goa_human.gaf.gz"

rule download_downloadGO:
	output:
		"source/go/go.json"
	shell:
		"cd transform/go && curl -o ../../source/go/go.json http://release.geneontology.org/2022-09-19/ontology/go.json"

rule go_gaf:
	input:
		"source/go/goa_human.gaf.gz",
		"transform/go/go_gaf.yaml"
	output:
		"output/go/go_gaf.dump.gaf.json.gz"
	shell:
		"sifter run transform/go/go_gaf.yaml"

rule goJson:
	input:
		"source/go/go.json",
		"schemas",
		"transform/go/go_json.yaml"
	output:
		"output/go/goJson.nodeOpen.before.json.gz",
		"output/go/goJson.nodeOpen.term.json.gz",
		"output/go/goJson.edgeOpen.edges.json.gz"
	shell:
		"sifter run transform/go/go_json.yaml"

rule goReduce:
	input:
		"output/go/go.edgeOpen.edges.json.gz",
		"transform/go/go_reduce.yaml"
	output:
		"output/go/goReduce.edgeReduce.edge_reduce.json.gz"
	shell:
		"sifter run transform/go/go_reduce.yaml"

rule download_gtex_downloadGeneTPM:
	output:
		"source/gtex/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz"
	shell:
		"cd transform/gtex && curl -o ../../source/gtex/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz"

rule download_gtex_downloadSampleAttributes:
	output:
		"source/gtex/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"
	shell:
		"cd transform/gtex && curl -o ../../source/gtex/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"

rule download_gtex_downloadSubjectPhenotypes:
	output:
		"source/gtex/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt"
	shell:
		"cd transform/gtex && curl -o ../../source/gtex/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt"

rule GTEX_Transcript_Expression:
	input:
		"source/gtex/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct.transpose.gz",
		"schema",
		"transform/gtex/transcript_transform.yaml"
	shell:
		"sifter run transform/gtex/transcript_transform.yaml"

rule mondo:
	input:
		"source/mondo/mondo.json",
		"transform/mondo/mondo.yaml"
	output:
		"output/mondo/mondo.extractEdges.edges.json.gz",
		"output/mondo/mondo.extractNodes.nodes.json.gz"
	shell:
		"sifter run transform/mondo/mondo.yaml"

rule downloads_msigdb:
	output:
		"source/msigdb/msigdb_v7.5.1.xml"
	shell:
		"cd transform/msigdb && curl -o ../../source/msigdb/msigdb_v7.5.1.xml https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.5.1/msigdb_v7.5.1.xml"

rule msigdb:
	input:
		"source/msigdb/msigdb_v7.5.1.xml",
		"schemas",
		"source/ensembl/gene2ensembl.translate.link.json.gz",
		"transform/msigdb/transform.yaml"
	output:
		"output/msigdb/msigdb.object.gene_set.json.gz",
		"output/msigdb/msigdb.transform.raw.json.gz",
		"output/msigdb/msigdb.transform.record.json.gz"
	shell:
		"sifter run transform/msigdb/transform.yaml"

rule pathwayCommonsDownload_pathwayCommonsDownload_pc_prep_download_drugbank:
	output:
		"source/pathway_commons/PathwayCommons12.drugbank.BIOPAX.owl.gz"
	shell:
		"cd transform/pathway_commons && curl -o ../../source/pathway_commons/PathwayCommons12.drugbank.BIOPAX.owl.gz https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.drugbank.BIOPAX.owl.gz"

rule pathwayCommonsDownload_pathwayCommonsDownload_pc_prep_download_hprd:
	output:
		"source/pathway_commons/PathwayCommons12.hprd.BIOPAX.owl.gz"
	shell:
		"cd transform/pathway_commons && curl -o ../../source/pathway_commons/PathwayCommons12.hprd.BIOPAX.owl.gz https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.hprd.BIOPAX.owl.gz"

rule pathwayCommonsDownload_pathwayCommonsDownload_pc_prep_process_drugbank:
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

rule pathwayCommonsDownload_pathwayCommonsDownload_pc_prep_process_reconx:
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

rule pathwayCommonsDownload_pathwayCommonsDownload_pc_prep_download_reconx:
	output:
		"source/pathway_commons/PathwayCommons12.reconx.BIOPAX.owl.gz"
	shell:
		"cd transform/pathway_commons && curl -o ../../source/pathway_commons/PathwayCommons12.reconx.BIOPAX.owl.gz https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.reconx.BIOPAX.owl.gz"

rule pathwayCommonsDownload_pathwayCommonsDownload_pc_prep_download_ctd:
	output:
		"source/pathway_commons/PathwayCommons12.ctd.BIOPAX.owl.gz"
	shell:
		"cd transform/pathway_commons && curl -o ../../source/pathway_commons/PathwayCommons12.ctd.BIOPAX.owl.gz https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.ctd.BIOPAX.owl.gz"

rule pathwayCommonsDownload_pathwayCommonsDownload_pc_prep_download_humancyc:
	output:
		"source/pathway_commons/PathwayCommons12.humancyc.BIOPAX.owl.gz"
	shell:
		"cd transform/pathway_commons && curl -o ../../source/pathway_commons/PathwayCommons12.humancyc.BIOPAX.owl.gz https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.humancyc.BIOPAX.owl.gz"

rule pathwayCommonsDownload_pathwayCommonsDownload_pc_prep_process_panther:
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

rule pathwayCommonsDownload_pathwayCommonsDownload_pc_prep_process_hprd:
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

rule pathwayCommonsDownload_pathwayCommonsDownload_pc_prep_download_dip:
	output:
		"source/pathway_commons/PathwayCommons12.dip.BIOPAX.owl.gz"
	shell:
		"cd transform/pathway_commons && curl -o ../../source/pathway_commons/PathwayCommons12.dip.BIOPAX.owl.gz https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.dip.BIOPAX.owl.gz"

rule pathwayCommonsDownload_pathwayCommonsDownload_pc_prep_download_pathbank:
	output:
		"source/pathway_commons/PathwayCommons12.pathbank.BIOPAX.owl.gz"
	shell:
		"cd transform/pathway_commons && curl -o ../../source/pathway_commons/PathwayCommons12.pathbank.BIOPAX.owl.gz https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.pathbank.BIOPAX.owl.gz"

rule pathwayCommonsDownload_pathwayCommonsDownload_pc_prep_download_netpath:
	output:
		"source/pathway_commons/PathwayCommons12.netpath.BIOPAX.owl.gz"
	shell:
		"cd transform/pathway_commons && curl -o ../../source/pathway_commons/PathwayCommons12.netpath.BIOPAX.owl.gz https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.netpath.BIOPAX.owl.gz"

rule pathwayCommonsDownload_pathwayCommonsDownload_pc_prep_process_dip:
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

rule pathwayCommonsDownload_pathwayCommonsDownload_pc_prep_process_psp:
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

rule pathwayCommonsDownload_pathwayCommonsDownload_pc_prep_download_corum:
	output:
		"source/pathway_commons/PathwayCommons12.corum.BIOPAX.owl.gz"
	shell:
		"cd transform/pathway_commons && curl -o ../../source/pathway_commons/PathwayCommons12.corum.BIOPAX.owl.gz https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.corum.BIOPAX.owl.gz"

rule pathwayCommonsDownload_pathwayCommonsDownload_pc_prep_download_msigdb:
	output:
		"source/pathway_commons/PathwayCommons12.msigdb.BIOPAX.owl.gz"
	shell:
		"cd transform/pathway_commons && curl -o ../../source/pathway_commons/PathwayCommons12.msigdb.BIOPAX.owl.gz https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.msigdb.BIOPAX.owl.gz"

rule pathwayCommonsDownload_pathwayCommonsDownload_pc_prep_download_biogrid:
	output:
		"source/pathway_commons/PathwayCommons12.biogrid.BIOPAX.owl.gz"
	shell:
		"cd transform/pathway_commons && curl -o ../../source/pathway_commons/PathwayCommons12.biogrid.BIOPAX.owl.gz https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.biogrid.BIOPAX.owl.gz"

rule pathwayCommonsDownload_pathwayCommonsDownload_pc_prep_download_pid:
	output:
		"source/pathway_commons/PathwayCommons12.pid.BIOPAX.owl.gz"
	shell:
		"cd transform/pathway_commons && curl -o ../../source/pathway_commons/PathwayCommons12.pid.BIOPAX.owl.gz https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.pid.BIOPAX.owl.gz"

rule pathwayCommonsDownload_pathwayCommonsDownload_pc_prep_download_mirtarbase:
	output:
		"source/pathway_commons/PathwayCommons12.mirtarbase.BIOPAX.owl.gz"
	shell:
		"cd transform/pathway_commons && curl -o ../../source/pathway_commons/PathwayCommons12.mirtarbase.BIOPAX.owl.gz https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.mirtarbase.BIOPAX.owl.gz"

rule pathwayCommonsDownload_pathwayCommonsDownload_pc_prep_download_reactome:
	output:
		"source/pathway_commons/PathwayCommons12.reactome.BIOPAX.owl.gz"
	shell:
		"cd transform/pathway_commons && curl -o ../../source/pathway_commons/PathwayCommons12.reactome.BIOPAX.owl.gz https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.reactome.BIOPAX.owl.gz"

rule pathwayCommonsDownload_pathwayCommonsDownload_pc_prep_process_reactome:
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

rule pathwayCommonsDownload_pathwayCommonsDownload_pc_prep_process_bind:
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

rule pathwayCommonsDownload_pathwayCommonsDownload_pc_prep_process_inoh:
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

rule pathwayCommonsDownload_pathwayCommonsDownload_pc_prep_download_bind:
	output:
		"source/pathway_commons/PathwayCommons12.bind.BIOPAX.owl.gz"
	shell:
		"cd transform/pathway_commons && curl -o ../../source/pathway_commons/PathwayCommons12.bind.BIOPAX.owl.gz https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.bind.BIOPAX.owl.gz"

rule pathwayCommonsDownload_pathwayCommonsDownload_pc_prep_download_panther:
	output:
		"source/pathway_commons/PathwayCommons12.panther.BIOPAX.owl.gz"
	shell:
		"cd transform/pathway_commons && curl -o ../../source/pathway_commons/PathwayCommons12.panther.BIOPAX.owl.gz https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.panther.BIOPAX.owl.gz"

rule pathwayCommonsDownload_pathwayCommonsDownload_pc_prep_download_innatedb:
	output:
		"source/pathway_commons/PathwayCommons12.innatedb.BIOPAX.owl.gz"
	shell:
		"cd transform/pathway_commons && curl -o ../../source/pathway_commons/PathwayCommons12.innatedb.BIOPAX.owl.gz https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.innatedb.BIOPAX.owl.gz"

rule pathwayCommonsDownload_pathwayCommonsDownload_pc_prep_process_pid:
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

rule pathwayCommonsDownload_pathwayCommonsDownload_pc_prep_process_corum:
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

rule pathwayCommonsDownload_pathwayCommonsDownload_pc_prep_process_humancyc:
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

rule pathwayCommonsDownload_pathwayCommonsDownload_pc_prep_process_innatedb:
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

rule pathwayCommonsDownload_pathwayCommonsDownload_pc_prep_download_psp:
	output:
		"source/pathway_commons/PathwayCommons12.psp.BIOPAX.owl.gz"
	shell:
		"cd transform/pathway_commons && curl -o ../../source/pathway_commons/PathwayCommons12.psp.BIOPAX.owl.gz https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.psp.BIOPAX.owl.gz"

rule pathwayCommonsDownload_jar:
	output:
		"transform/pathway_commons/target/pc-extract-1.0-jar-with-dependencies.jar"
	shell:
		"cd transform/pathway_commons && docker run  --rm -u `id -u` -v `pwd`:`pwd` -v `pwd`/maven-repo:/root/.m2 -w `pwd` -ti maven:3-openjdk-11 mvn package"

rule pathwayCommonsDownload_pathwayCommonsDownload_pc_prep_process_mirtarbase:
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

rule pathwayCommonsDownload_pathwayCommonsDownload_pc_prep_download_inoh:
	output:
		"source/pathway_commons/PathwayCommons12.inoh.BIOPAX.owl.gz"
	shell:
		"cd transform/pathway_commons && curl -o ../../source/pathway_commons/PathwayCommons12.inoh.BIOPAX.owl.gz https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.inoh.BIOPAX.owl.gz"

rule pathwayCommonsDownload_pathwayCommonsDownload_pc_prep_process_kegg:
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

rule pathwayCommonsDownload_pathwayCommonsDownload_pc_prep_process_biogrid:
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

rule pathwayCommonsDownload_pathwayCommonsDownload_pc_prep_process_msigdb:
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

rule pathwayCommonsDownload_pathwayCommonsDownload_pc_prep_download_kegg:
	output:
		"source/pathway_commons/PathwayCommons12.kegg.BIOPAX.owl.gz"
	shell:
		"cd transform/pathway_commons && curl -o ../../source/pathway_commons/PathwayCommons12.kegg.BIOPAX.owl.gz https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.kegg.BIOPAX.owl.gz"

rule pathwayCommonsDownload_pathwayCommonsDownload_pc_prep_process_pathbank:
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

rule pathwayCommonsDownload_pathwayCommonsDownload_pc_prep_process_ctd:
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

rule pathwayCommonsDownload_pathwayCommonsDownload_pc_prep_process_netpath:
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

rule pathway_commons:
	input:
		"source/pathway_commons",
		"source/pathway_commons",
		"schema",
		"transform/pathway_commons/transform.yaml"
	output:
		"output/pathway_commons/pathway_commons.complexBundle.complex.json.gz",
		"output/pathway_commons/pathway_commons.interactionMap.interaction.json.gz"
	shell:
		"sifter run transform/pathway_commons/transform.yaml"

rule pharmacodb:
	input:
		"source/pharmacodb/pharmacodb-1.1.1.sql",
		"transform/pharmacodb/database.yaml"
	output:
		"source/pharmacodb/tables/source_tissue_names.tsv.gz",
		"source/pharmacodb/tables/cell_tissues.tsv.gz",
		"source/pharmacodb/tables/dose_responses.tsv.gz",
		"source/pharmacodb/tables/sources.tsv.gz",
		"source/pharmacodb/tables/cellosaurus.tsv.gz",
		"source/pharmacodb/tables/datasets.tsv.gz",
		"source/pharmacodb/tables/source_drug_names.tsv.gz",
		"source/pharmacodb/tables/cells.tsv.gz",
		"source/pharmacodb/tables/dataset_cells.tsv.gz",
		"source/pharmacodb/tables/drugs.tsv.gz",
		"source/pharmacodb/tables/source_statistics.tsv.gz",
		"source/pharmacodb/tables/source_cell_names.tsv.gz",
		"source/pharmacodb/tables/drug_annots.tsv.gz",
		"source/pharmacodb/tables/tissues.tsv.gz",
		"source/pharmacodb/tables/experiments.tsv.gz",
		"source/pharmacodb/tables/profiles.tsv.gz"
	shell:
		"sifter run transform/pharmacodb/database.yaml"

rule pharmacodbDownload_download:
	output:
		"source/pharmacodb/pharmacodb-1.1.1.sql"
	shell:
		"cd transform/pharmacodb && curl -o ../../source/pharmacodb/pharmacodb-1.1.1.sql https://zenodo.org/record/1143645/files/pharmacodb-1.1.1.sql"

rule pharmacodb_profiles:
	input:
		"source/pharmacodb/tables/profiles.tsv.gz",
		"schema",
		"schema",
		"schema",
		"schema",
		"schema",
		"source/pharmacodb/tables/dr_reduce.curveReduce.dose_response_curve.json.gz",
		"source/pharmacodb/tables/experiments.tsv.gz",
		"source/pharmacodb/tables/drugs.tsv.gz",
		"source/pharmacodb/tables/datasets.tsv.gz",
		"source/pharmacodb/tables/cells.tsv.gz",
		"transform/pharmacodb/profile.yaml"
	output:
		"output/pharmacodb/pharmacodb_profiles.cellAliquot.aliquot.json.gz",
		"output/pharmacodb/pharmacodb_profiles.cellCase.case.json.gz",
		"output/pharmacodb/pharmacodb_profiles.cellDistinct.checkpoint.json.gz",
		"output/pharmacodb/pharmacodb_profiles.cellProject.project.json.gz",
		"output/pharmacodb/pharmacodb_profiles.cellSample.sample.json.gz",
		"output/pharmacodb/pharmacodb_profiles.drObject.drug_response.json.gz"
	shell:
		"sifter run transform/pharmacodb/profile.yaml"

rule dr_reduce:
	input:
		"source/pharmacodb/tables/dose_responses.tsv.gz",
		"transform/pharmacodb/reduce_dr_data.yaml"
	output:
		"source/pharmacodb/tables/dr_reduce.curveReduce.dose_response_curve.json.gz"
	shell:
		"sifter run transform/pharmacodb/reduce_dr_data.yaml"

rule pubmed:
	input:
		"source/pubmed/baseline",
		"schema",
		"transform/pubmed/transform.yaml"
	output:
		"output/pubmed/pubmed.transform.publication.json.gz"
	shell:
		"sifter run transform/pubmed/transform.yaml"



