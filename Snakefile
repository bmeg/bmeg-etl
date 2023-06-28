
scattergather:
	compoundDist_compoundDist=500

rule all:
	input:
		"graph/graphbuild_pdb.start-graph.vertex.json.gz",
		"graph/graphbuild_gdc.projectObject-graph.vertex.json.gz",
		"graph/graphbuild_mondo.extract-graph.vertex.json.gz",
		"graph/graphbuild_gdc-mafs.somaticCallsets-graph.edge.json.gz",
		"graph/graphbuild_GDSC_Transform.transform-graph.vertex.json.gz",
		"graph/graphbuild_chemblDrugMechanismExtract.build-graph.edge.json.gz",
		"graph/graphbuild_gdc.caseObject-graph.vertex.json.gz",
		"graph/graphbuild_depmap-cases.aliquots-graph.edge.json.gz",
		"graph/graphbuild_depmap-mafs.callsets-graph.edge.json.gz",
		"tables/pharmacodb/dr_reduce.curveReduce.dose_response_curve.json.gz",
		"graph/graphbuild_go.transform-graph.edge.json.gz",
		"graph/graphbuild_pubmed.transform-graph.edge.json.gz",
		"graph/graphbuild_ensembl_gtf.transcripts-graph.vertex.json.gz",
		"graph/graphbuild_pharmacodb_profiles_gx.drugResponse-graph.edge.json.gz",
		"output/pharmacodb/pharmacodb_profiles_gx.transform.debug.json.gz",
		"graph/graphbuild_GDSC_Transform.transform-graph.edge.json.gz",
		"graph/graphbuild_ensembl_gtf.exons-graph.vertex.json.gz",
		"graph/graphbuild_cellosaurus.transform-graph.vertex.json.gz",
		"graph/graphbuild_cellosaurus.transform-graph.edge.json.gz",
		"graph/graphbuild_g2p.main-graph.edge.json.gz",
		"graph/graphbuild_GTEX_Transcript_Expression.update-graph.vertex.json.gz",
		"graph/graphbuild_depmap-mafs.variants-graph.edge.json.gz",
		"graph/graphbuild_gdc.aliquotObject-graph.vertex.json.gz",
		"graph/graphbuild_go.transform-graph.vertex.json.gz",
		"graph/graphbuild_pharmacodb_profiles_gx.sample-graph.vertex.json.gz",
		"graph/graphbuild_GTEX_Gene_Expression.gctProcess-graph.vertex.json.gz",
		"graph/graphbuild_GTEX_Transcript_Expression.update-graph.edge.json.gz",
		"graph/graphbuild_rnaseq.rna-graph.vertex.json.gz",
		"graph/graphbuild_pharmacodb_profiles_gx.sample-graph.edge.json.gz",
		"graph/graphbuild_gdc.sampleObject-graph.vertex.json.gz",
		"graph/graphbuild_annotatedAllele.alleleEffect-graph.vertex.json.gz",
		"graph/graphbuild_gdc-mafs.somaticCallsets-graph.vertex.json.gz",
		"graph/graphbuild_pathway_commons.interactionMap-graph.vertex.json.gz",
		"graph/graphbuild_pubmed.transform-graph.vertex.json.gz",
		"graph/graphbuild_depmap-expression.values-graph.edge.json.gz",
		"graph/graphbuild_depmap-mafs.callsets-graph.vertex.json.gz",
		"graph/graphbuild_depmap-mafs.variants-graph.vertex.json.gz",
		"output/go/go_gaf.dump.gaf.json.gz",
		"graph/graphbuild_annotatedAllele.allele-graph.edge.json.gz",
		"graph/graphbuild_chemblTransform.records-graph.edge.json.gz",
		"graph/graphbuild_uniprot_trembl.start-graph.edge.json.gz",
		"graph/graphbuild_uniprot_sprot.start-graph.edge.json.gz",
		"normalize/chembl/compounds.dist.compoundDistant.json.gz",
		"graph/graphbuild_pharmacodb_profiles_gx.drugResponse-graph.vertex.json.gz",
		"graph/graphbuild_prism_transform.secondary-graph.edge.json.gz",
		"graph/graphbuild_depmap-cases.aliquots-graph.vertex.json.gz",
		"graph/graphbuild_depmap-expression.values-graph.vertex.json.gz",
		"graph/graphbuild_GDSC_VCF_Transform.variants-graph.vertex.json.gz",
		"graph/graphbuild_msigdb.transform-graph.edge.json.gz",
		"graph/graphbuild_prism_transform.primary-graph.edge.json.gz",
		"graph/graphbuild_annotatedAllele.alleleEffect-graph.edge.json.gz",
		"graph/graphbuild_gdc.aliquotObject-graph.edge.json.gz",
		"graph/graphbuild_gdc.caseObject-graph.edge.json.gz",
		"graph/graphbuild_ensembl_gtf.genes-graph.vertex.json.gz",
		"graph/graphbuild_pathway_commons.interactionMap-graph.edge.json.gz",
		"output/gdc/gdc.aliquotAlias.table.json.gz",
		"graph/graphbuild_GDSC_VCF_Transform.callset-graph.edge.json.gz",
		"graph/graphbuild_GDSC_rnaseq_Transform.start-graph.vertex.json.gz",
		"graph/graphbuild_bindingdbTsv.start-graph.edge.json.gz",
		"graph/graphbuild_ensembl_gtf.exons-graph.edge.json.gz",
		"graph/graphbuild_prism_transform.primary-graph.vertex.json.gz",
		"graph/graphbuild_uniprot_trembl.start-graph.vertex.json.gz",
		"graph/graphbuild_chemblDrugMechanismExtract.build-graph.vertex.json.gz",
		"graph/graphbuild_pdb.start-graph.edge.json.gz",
		"output/pathway_commons/pathway_commons.complexBundle.complex.json.gz",
		"output-normalize/allele.vcf",
		"graph/graphbuild_annotatedAllele.allele-graph.vertex.json.gz",
		"graph/graphbuild_pharmacodb_profiles_gx.aliquot-graph.vertex.json.gz",
		"graph/graphbuild_chemblTransform.records-graph.vertex.json.gz",
		"graph/graphbuild_ensembl_gtf.transcripts-graph.edge.json.gz",
		"graph/graphbuild_gdc-mafs.scan-graph.edge.json.gz",
		"graph/graphbuild_msigdb.transform-graph.vertex.json.gz",
		"graph/graphbuild_depmap-cases.depMapCases-graph.vertex.json.gz",
		"graph/graphbuild_GTEX_Gene_Expression.gctProcess-graph.edge.json.gz",
		"graph/graphbuild_gdc.sampleObject-graph.edge.json.gz",
		"graph/graphbuild_prism_transform.secondary-graph.vertex.json.gz",
		"graph/graphbuild_bindingdbTsv.start-graph.vertex.json.gz",
		"graph/graphbuild_GDSC_rnaseq_Transform.aliquot-graph.vertex.json.gz",
		"graph/graphbuild_ensembl_gtf.genes-graph.edge.json.gz",
		"graph/graphbuild_gdc-mafs.scan-graph.vertex.json.gz",
		"graph/graphbuild_GDSC_rnaseq_Transform.start-graph.edge.json.gz",
		"graph/graphbuild_GDSC_Transform.aliquot-graph.vertex.json.gz",
		"graph/graphbuild_GDSC_VCF_Transform.callset-graph.vertex.json.gz",
		"graph/graphbuild_GDSC_VCF_Transform.variants-graph.edge.json.gz",
		"graph/graphbuild_depmap-cases.depMapCases-graph.edge.json.gz",
		"graph/graphbuild_pharmacodb_profiles_gx.aliquot-graph.edge.json.gz",
		"source/g2p/tables/g2p_source.bed",
		"graph/graphbuild_g2p.main-graph.vertex.json.gz",
		"graph/graphbuild_gdc.projectObject-graph.edge.json.gz",
		"graph/graphbuild_mondo.extract-graph.edge.json.gz",
		"graph/graphbuild_uniprot_sprot.start-graph.vertex.json.gz",
		"output/gdsc/GDSC_VCF_Transform.namefix.source.json.gz",
		"graph/graphbuild_depmap-cases.sampleObjects-graph.edge.json.gz",
		"graph/graphbuild_GDSC_rnaseq_Transform.aliquot-graph.edge.json.gz",
		"graph/graphbuild_depmap-cases.sampleObjects-graph.vertex.json.gz",
		"graph/graphbuild_rnaseq.rna-graph.edge.json.gz",
		"graph/graphbuild_GDSC_Transform.aliquot-graph.edge.json.gz"

rule bindingdbTsv:
	input:
		"source/bindingdb/BindingDB_All.tsv",
		"schema",
		"transform/bindingdb/transform_tsv.yaml"
	output:
		"output/bindingdb/bindingdbTsv.start.protein_compound_association.json.gz"
	shell:
		"sifter run transform/bindingdb/transform_tsv.yaml"

rule cellosarusSynonyms:
	input:
		"source/cellosaurus/cellosaurus.obo",
		"transform/cellosaurus/alias_table.yaml"
	output:
		"tables/cellosarusSynonyms.caseTable.ach2cellosaurus.json.gz"
	shell:
		"sifter run transform/cellosaurus/alias_table.yaml"

rule cellosaurus:
	input:
		"source/cellosaurus/cellosaurus.obo",
		"tables/ncit2mondo.ncit_extract.mapping.json.gz",
		"schema",
		"transform/cellosaurus/transform.yaml"
	output:
		"output/cellosaurus/cellosaurus.transform.cases.json.gz"
	shell:
		"sifter run transform/cellosaurus/transform.yaml"

rule compoundDist_compoundDist_compoundDist_gather:
	input:
		input= gather.compoundDist_compoundDist("transform/chembl/shards/{scatteritem}")
	output:
		output= "normalize/chembl/compounds.dist.compoundDistant.json.gz"
	shell:
		"cat {input} > {output}"

rule compoundDist_compoundDist_compoundDist_scatter:
	input:
		input= "output/chembl/chemblTransform.records.compound.json.gz"
	output:
		output= "transform/chembl/shards/{shard}-of-{total}"
	resources:
		mem_mb=30000
	shell:
		"/home/groups/EllrottLab/bmeg-etl/util/compound_distance.py -n {threads} -i {input} -s {wildcards.shard} -t {wildcards.total} -o {output}"

rule chemblDrugMechanismExtract:
	input:
		"source/chembl/chembl_31/chembl_31_sqlite/chembl_31.db",
		"schema",
		"transform/chembl/drug_mechanism.yaml"
	output:
		"output/chembl/chemblDrugMechanismExtract.build.protein_compound_association.json.gz"
	shell:
		"sifter run transform/chembl/drug_mechanism.yaml"

rule chemblSynonyms:
	input:
		"source/chembl/chembl_31/chembl_31_sqlite/chembl_31.db",
		"transform/chembl/synonyms.yaml"
	output:
		"tables/chemblSynonyms.buildTable.synonyms.json.gz",
		"tables/chemblSynonyms.longTable.synonyms.json.gz"
	shell:
		"sifter run transform/chembl/synonyms.yaml"

rule chemblTransform:
	input:
		"source/chembl/chembl_31/chembl_31_sqlite/chembl_31.db",
		"tables/chemblSynonyms.buildTable.synonyms.json.gz",
		"schema",
		"transform/chembl/transform.yaml"
	output:
		"output/chembl/chemblTransform.records.compound.json.gz"
	shell:
		"sifter run transform/chembl/transform.yaml"

rule depmap_cases:
	input:
		"source/depmap/Model.csv",
		"schema",
		"tables/cellosarusSynonyms.caseTable.ach2cellosaurus.json.gz",
		"schema",
		"schema",
		"transform/depmap/cases.yaml"
	output:
		"output/depmap/depmap-cases.aliquots.aliquot.json.gz",
		"output/depmap/depmap-cases.depMapCases.case.json.gz",
		"output/depmap/depmap-cases.sampleObjects.sample.json.gz"
	shell:
		"sifter run transform/depmap/cases.yaml"

rule depmap_expression:
	input:
		"source/depmap/OmicsExpressionProteinCodingGenesTPMLogp1.csv",
		"tables/gene2ensembl.translate.link.json.gz",
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
		"schema",
		"transform/depmap/mutations.yaml"
	output:
		"output/depmap/depmap-mafs.allele.allele.json.gz",
		"output/depmap/depmap-mafs.callsets.callset.json.gz",
		"output/depmap/depmap-mafs.variants.variants.json.gz"
	shell:
		"sifter run transform/depmap/mutations.yaml"

rule ensembl_gtf:
	input:
		"source/ensembl/Homo_sapiens.GRCh38.108.chr_patch_hapl_scaff.gtf.gz",
		"schema",
		"schema",
		"schema",
		"transform/ensembl/ensembl_transform.yaml"
	output:
		"output/ensembl/ensembl_gtf.genes.gene.json.gz",
		"output/ensembl/ensembl_gtf.transcripts.transcript.json.gz",
		"output/ensembl/ensembl_gtf.exons.exon.json.gz"
	shell:
		"sifter run transform/ensembl/ensembl_transform.yaml"

rule gene2ensembl:
	input:
		"source/ensembl/gene2ensembl.gz",
		"transform/ensembl/gene2ensembl.yaml"
	output:
		"tables/gene2ensembl.translate.link.json.gz"
	shell:
		"sifter run transform/ensembl/gene2ensembl.yaml"

rule hugo2ensembl:
	input:
		"source/ensembl/Homo_sapiens.GRCh38.108.chr_patch_hapl_scaff.gtf.gz",
		"transform/ensembl/hugo2ensembl.yaml"
	output:
		"tables/hugo2ensembl.tsv"
	shell:
		"sifter run transform/ensembl/hugo2ensembl.yaml"

rule g2p_bedfile:
	input:
		"source/g2p",
		"transform/g2p/bed_file.yaml"
	output:
		"source/g2p/tables/g2p_source.bed"
	shell:
		"sifter run transform/g2p/bed_file.yaml"

rule g2p:
	input:
		"source/g2p",
		"schema",
		"tables/hugo2ensembl.tsv",
		"source/g2p/tables/hglft_genome_2749d_8437f0.bed",
		"schema",
		"transform/g2p/transform.yaml"
	output:
		"output/g2p/g2p.alleles.allele.json.gz",
		"output/g2p/g2p.main.assocation.json.gz"
	shell:
		"sifter run transform/g2p/transform.yaml"

rule gdc_mafs:
	input:
		"source/gdc/open-maf",
		"schema",
		"schema",
		"schema",
		"transform/gdc/maf-files.yaml"
	output:
		"output/gdc/open-maf/gdc-mafs.allele.allele.json.gz",
		"output/gdc/open-maf/gdc-mafs.scan.variant.json.gz",
		"output/gdc/open-maf/gdc-mafs.somaticCallsets.callset.json.gz"
	shell:
		"sifter run transform/gdc/maf-files.yaml"

rule rnaseq:
	input:
		"source/gdc/files.json",
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
		"schema",
		"transform/gdc/transform.yaml"
	output:
		"output/gdc/gdc.caseObject.case.json.gz",
		"output/gdc/gdc.projectObject.project.json.gz",
		"output/gdc/gdc.sampleObject.sample.json.gz",
		"output/gdc/gdc.aliquotAlias.table.json.gz",
		"output/gdc/gdc.aliquotObject.aliquot.json.gz"
	shell:
		"sifter run transform/gdc/transform.yaml"

rule cosmic2ach:
	input:
		"source/gdsc/model_list_20191104.csv",
		"transform/gdsc/cosmic_to_ach.yaml"
	output:
		"tables/cosmic2ach.translate.link.json.gz"
	shell:
		"sifter run transform/gdsc/cosmic_to_ach.yaml"

rule GDSC_Transform:
	input:
		"source/gdsc/Cell_line_RMA_proc_basalExp.txt",
		"tables/cosmic2ach.translate.link.json.gz",
		"schema",
		"tables/hugo2ensembl.tsv",
		"schema",
		"transform/gdsc/rna_transform.yaml"
	output:
		"output/gdsc/GDSC_Transform.aliquot.aliquot.json.gz",
		"output/gdsc/GDSC_Transform.transform.geneExpression.json.gz"
	shell:
		"sifter run transform/gdsc/rna_transform.yaml"

rule GDSC_rnaseq_Transform:
	input:
		"source/gdsc/rnaseq_sanger_20210316.csv",
		"source/gdsc/model_list_20230110.csv",
		"schema",
		"source/gdsc/gene_identifiers_20191101.csv",
		"schema",
		"transform/gdsc/rnaseq_transform.yaml"
	output:
		"output/gdsc/GDSC_rnaseq_Transform.aliquot.aliquot.json.gz",
		"output/gdsc/GDSC_rnaseq_Transform.start.geneExpression.json.gz"
	shell:
		"sifter run transform/gdsc/rnaseq_transform.yaml"

rule GDSC_VCF_Transform:
	input:
		"source/gdsc/mutations_wes_vcf_20221010",
		"schema",
		"schema",
		"source/gdsc/model_list_20230110.csv",
		"schema",
		"transform/gdsc/vcf_transform.yaml"
	output:
		"output/gdsc/GDSC_VCF_Transform.allele.allele.json.gz",
		"output/gdsc/GDSC_VCF_Transform.callset.somatic_callset.json.gz",
		"output/gdsc/GDSC_VCF_Transform.namefix.source.json.gz",
		"output/gdsc/GDSC_VCF_Transform.variants.somatic_variant.json.gz"
	shell:
		"sifter run transform/gdsc/vcf_transform.yaml"

rule go_gaf:
	input:
		"source/go/goa_human.gaf.gz",
		"transform/go/go_gaf.yaml"
	output:
		"output/go/go_gaf.dump.gaf.json.gz"
	shell:
		"sifter run transform/go/go_gaf.yaml"

rule go:
	input:
		"source/go/go.obo",
		"schema",
		"transform/go/go_obo.yaml"
	output:
		"output/go/go.transform.term.json.gz"
	shell:
		"sifter run transform/go/go_obo.yaml"

rule GTEX_Gene_Expression:
	input:
		"source/gtex/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz",
		"schema",
		"transform/gtex/gene_transform.yaml"
	output:
		"output/gtex/GTEX_Gene_Expression.gctProcess.gene_expression.json.gz"
	resources:
		mem_mb=50000
	shell:
		"sifter run transform/gtex/gene_transform.yaml"

rule GTEX_Transcript_Expression:
	input:
		"source/gtex/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct.gz",
		"schema",
		"transform/gtex/transcript_transform.yaml"
	output:
		"output/gtex/GTEX_Transcript_Expression.update.transcript_expression.json.gz"
	resources:
		mem_mb=80000
	shell:
		"sifter run transform/gtex/transcript_transform.yaml"

rule interpro:
	input:
		"source/interpro/interpro.xml.gz",
		"transform/interpro/transform.yaml"
	shell:
		"sifter run transform/interpro/transform.yaml"

rule mondo:
	input:
		"source/mondo/mondo.json",
		"schema",
		"transform/mondo/mondo.yaml"
	output:
		"output/mondo/mondo.extract.phenotype.json.gz"
	shell:
		"sifter run transform/mondo/mondo.yaml"

rule ncit2mondo:
	input:
		"source/mondo/mondo.json",
		"transform/mondo/ncit2mondo_table.yaml"
	output:
		"tables/ncit2mondo.ncit_extract.mapping.json.gz"
	shell:
		"sifter run transform/mondo/ncit2mondo_table.yaml"

rule msigdb:
	input:
		"source/msigdb/msigdb_v7.5.1.xml",
		"tables/gene2ensembl.translate.link.json.gz",
		"schema",
		"transform/msigdb/transform.yaml"
	output:
		"output/msigdb/msigdb.transform.gene_set.json.gz"
	shell:
		"sifter run transform/msigdb/transform.yaml"

rule pathway_commons:
	input:
		"source/pathway_commons",
		"source/pathway_commons",
		"schema",
		"transform/pathway_commons/transform.yaml"
	output:
		"output/pathway_commons/pathway_commons.complexBundle.complex.json.gz",
		"output/pathway_commons/pathway_commons.interactionMap.interaction.json.gz"
	resources:
		mem_mb=20000
	shell:
		"sifter run transform/pathway_commons/transform.yaml"

rule pdb:
	input:
		"source/pdb/entries.idx",
		"schema",
		"transform/pdb/transform.yaml"
	output:
		"output/pdb/pdb.start.protein_structure.json.gz"
	shell:
		"sifter run transform/pdb/transform.yaml"

rule pharmacoTableGenerate_pharmacodb:
	input:
		script= "transform/pharmacodb/drug_response_tsv.R",
		source= "source/pharmacodb"
	output:
		profile= "tables/pharmacodb_profiles.tsv"
	shell:
		"{input.script} -o {output.profile} -s {input.source}"

rule pharmacodb_profiles_gx:
	input:
		"tables/pharmacodb_profiles.tsv",
		"tables/chemblSynonyms.longTable.synonyms.json.gz",
		"tables/chemblSynonyms.longTable.synonyms.json.gz",
		"tables/chemblSynonyms.longTable.synonyms.json.gz",
		"schema",
		"schema",
		"schema",
		"transform/pharmacodb/profile.yaml"
	output:
		"output/pharmacodb/pharmacodb_profiles_gx.transform.debug.json.gz",
		"output/pharmacodb/pharmacodb_profiles_gx.aliquot.aliquot.json.gz",
		"output/pharmacodb/pharmacodb_profiles_gx.drugResponse.drug_response.json.gz",
		"output/pharmacodb/pharmacodb_profiles_gx.sample.sample.json.gz"
	shell:
		"sifter run transform/pharmacodb/profile.yaml"

rule dr_reduce:
	input:
		"tables/pharmacodb/dose_responses.tsv.gz",
		"transform/pharmacodb/reduce_dr_data.yaml"
	output:
		"tables/pharmacodb/dr_reduce.curveReduce.dose_response_curve.json.gz"
	shell:
		"sifter run transform/pharmacodb/reduce_dr_data.yaml"

rule prism_transform:
	input:
		"source/prism/primary-screen-replicate-collapsed-logfold-change.csv",
		"source/prism/secondary-screen-replicate-collapsed-logfold-change.csv",
		"source/prism/primary-screen-replicate-collapsed-treatment-info.csv",
		"schema",
		"source/prism/secondary-screen-replicate-collapsed-treatment-info.csv",
		"schema",
		"transform/prism/transform.yaml"
	output:
		"output/prism/prism_transform.primary.drug_response.json.gz",
		"output/prism/prism_transform.secondary.drug_response.json.gz"
	shell:
		"sifter run transform/prism/transform.yaml"

rule pubmed:
	input:
		"source/pubmed/baseline",
		"schema",
		"transform/pubmed/transform.yaml"
	output:
		"output/pubmed/pubmed.transform.publication.json.gz"
	resources:
		mem_mb=20000
	shell:
		"sifter run transform/pubmed/transform.yaml"

rule uniprot_sprot:
	input:
		"source/uniprot/uniprot_sprot_human.xml.gz",
		"schema",
		"transform/uniprot/transform_sprot.yaml"
	output:
		"output/uniprot/uniprot_sprot.start.protein.json.gz"
	shell:
		"sifter run transform/uniprot/transform_sprot.yaml"

rule uniprot_trembl:
	input:
		"source/uniprot/uniprot_trembl_human.xml.gz",
		"schema",
		"transform/uniprot/transform_trembl.yaml"
	output:
		"output/uniprot/uniprot_trembl.start.protein.json.gz"
	shell:
		"sifter run transform/uniprot/transform_trembl.yaml"

rule collect_Allele:
	input:
		"output/depmap/depmap-mafs.allele.allele.json.gz",
		"output/g2p/g2p.alleles.allele.json.gz",
		"output/gdc/open-maf/gdc-mafs.allele.allele.json.gz",
		"output/gdsc/GDSC_VCF_Transform.allele.allele.json.gz"
	output:
		"output-normalize/allele.merge.json.gz"
	shell:
		"lathe class-concat Allele transform -o output-normalize/allele.merge.json.gz"

rule annotatedAllele:
	input:
		"output-normalize/allele.annotated.vcf",
		"schema",
		"schema",
		"normalize/allele/annotate_transform.yaml"
	output:
		"output/allele/annotatedAllele.allele.allele.json.gz",
		"output/allele/annotatedAllele.alleleEffect.alleleEffect.json.gz"
	shell:
		"sifter run normalize/allele/annotate_transform.yaml"

rule annotated_VCF_transform:
	input:
		"output-normalize/allele.merge.json.gz",
		"normalize/allele/vcf_transform.yaml"
	output:
		"output-normalize/allele.vcf"
	shell:
		"sifter run normalize/allele/vcf_transform.yaml"

rule graphbuild_GDSC_Transform:
	input:
		"output/gdsc/GDSC_Transform.aliquot.aliquot.json.gz",
		"output/gdsc/GDSC_Transform.transform.geneExpression.json.gz",
		"schema",
		"schema",
		"graph-build/graphbuild_GDSC_Transform.yaml"
	output:
		"graph/graphbuild_GDSC_Transform.aliquot-graph.vertex.json.gz",
		"graph/graphbuild_GDSC_Transform.aliquot-graph.edge.json.gz",
		"graph/graphbuild_GDSC_Transform.transform-graph.vertex.json.gz",
		"graph/graphbuild_GDSC_Transform.transform-graph.edge.json.gz"
	shell:
		"sifter run graph-build/graphbuild_GDSC_Transform.yaml"

rule graphbuild_GDSC_VCF_Transform:
	input:
		"output/gdsc/GDSC_VCF_Transform.callset.somatic_callset.json.gz",
		"output/gdsc/GDSC_VCF_Transform.variants.somatic_variant.json.gz",
		"schema",
		"schema",
		"graph-build/graphbuild_GDSC_VCF_Transform.yaml"
	output:
		"graph/graphbuild_GDSC_VCF_Transform.callset-graph.vertex.json.gz",
		"graph/graphbuild_GDSC_VCF_Transform.callset-graph.edge.json.gz",
		"graph/graphbuild_GDSC_VCF_Transform.variants-graph.vertex.json.gz",
		"graph/graphbuild_GDSC_VCF_Transform.variants-graph.edge.json.gz"
	shell:
		"sifter run graph-build/graphbuild_GDSC_VCF_Transform.yaml"

rule graphbuild_GDSC_rnaseq_Transform:
	input:
		"output/gdsc/GDSC_rnaseq_Transform.aliquot.aliquot.json.gz",
		"output/gdsc/GDSC_rnaseq_Transform.start.geneExpression.json.gz",
		"schema",
		"schema",
		"graph-build/graphbuild_GDSC_rnaseq_Transform.yaml"
	output:
		"graph/graphbuild_GDSC_rnaseq_Transform.aliquot-graph.vertex.json.gz",
		"graph/graphbuild_GDSC_rnaseq_Transform.aliquot-graph.edge.json.gz",
		"graph/graphbuild_GDSC_rnaseq_Transform.start-graph.vertex.json.gz",
		"graph/graphbuild_GDSC_rnaseq_Transform.start-graph.edge.json.gz"
	shell:
		"sifter run graph-build/graphbuild_GDSC_rnaseq_Transform.yaml"

rule graphbuild_GTEX_Gene_Expression:
	input:
		"output/gtex/GTEX_Gene_Expression.gctProcess.gene_expression.json.gz",
		"schema",
		"graph-build/graphbuild_GTEX_Gene_Expression.yaml"
	output:
		"graph/graphbuild_GTEX_Gene_Expression.gctProcess-graph.vertex.json.gz",
		"graph/graphbuild_GTEX_Gene_Expression.gctProcess-graph.edge.json.gz"
	shell:
		"sifter run graph-build/graphbuild_GTEX_Gene_Expression.yaml"

rule graphbuild_GTEX_Transcript_Expression:
	input:
		"output/gtex/GTEX_Transcript_Expression.update.transcript_expression.json.gz",
		"schema",
		"graph-build/graphbuild_GTEX_Transcript_Expression.yaml"
	output:
		"graph/graphbuild_GTEX_Transcript_Expression.update-graph.edge.json.gz",
		"graph/graphbuild_GTEX_Transcript_Expression.update-graph.vertex.json.gz"
	shell:
		"sifter run graph-build/graphbuild_GTEX_Transcript_Expression.yaml"

rule graphbuild_annotatedAllele:
	input:
		"output/allele/annotatedAllele.allele.allele.json.gz",
		"output/allele/annotatedAllele.alleleEffect.alleleEffect.json.gz",
		"schema",
		"schema",
		"graph-build/graphbuild_annotatedAllele.yaml"
	output:
		"graph/graphbuild_annotatedAllele.allele-graph.vertex.json.gz",
		"graph/graphbuild_annotatedAllele.allele-graph.edge.json.gz",
		"graph/graphbuild_annotatedAllele.alleleEffect-graph.vertex.json.gz",
		"graph/graphbuild_annotatedAllele.alleleEffect-graph.edge.json.gz"
	shell:
		"sifter run graph-build/graphbuild_annotatedAllele.yaml"

rule graphbuild_bindingdbTsv:
	input:
		"output/bindingdb/bindingdbTsv.start.protein_compound_association.json.gz",
		"schema",
		"graph-build/graphbuild_bindingdbTsv.yaml"
	output:
		"graph/graphbuild_bindingdbTsv.start-graph.vertex.json.gz",
		"graph/graphbuild_bindingdbTsv.start-graph.edge.json.gz"
	shell:
		"sifter run graph-build/graphbuild_bindingdbTsv.yaml"

rule graphbuild_cellosaurus:
	input:
		"output/cellosaurus/cellosaurus.transform.cases.json.gz",
		"schema",
		"graph-build/graphbuild_cellosaurus.yaml"
	output:
		"graph/graphbuild_cellosaurus.transform-graph.vertex.json.gz",
		"graph/graphbuild_cellosaurus.transform-graph.edge.json.gz"
	shell:
		"sifter run graph-build/graphbuild_cellosaurus.yaml"

rule graphbuild_chemblDrugMechanismExtract:
	input:
		"output/chembl/chemblDrugMechanismExtract.build.protein_compound_association.json.gz",
		"schema",
		"graph-build/graphbuild_chemblDrugMechanismExtract.yaml"
	output:
		"graph/graphbuild_chemblDrugMechanismExtract.build-graph.vertex.json.gz",
		"graph/graphbuild_chemblDrugMechanismExtract.build-graph.edge.json.gz"
	shell:
		"sifter run graph-build/graphbuild_chemblDrugMechanismExtract.yaml"

rule graphbuild_chemblTransform:
	input:
		"output/chembl/chemblTransform.records.compound.json.gz",
		"schema",
		"graph-build/graphbuild_chemblTransform.yaml"
	output:
		"graph/graphbuild_chemblTransform.records-graph.edge.json.gz",
		"graph/graphbuild_chemblTransform.records-graph.vertex.json.gz"
	shell:
		"sifter run graph-build/graphbuild_chemblTransform.yaml"

rule graphbuild_depmap_cases:
	input:
		"output/depmap/depmap-cases.sampleObjects.sample.json.gz",
		"output/depmap/depmap-cases.aliquots.aliquot.json.gz",
		"output/depmap/depmap-cases.depMapCases.case.json.gz",
		"schema",
		"schema",
		"schema",
		"graph-build/graphbuild_depmap-cases.yaml"
	output:
		"graph/graphbuild_depmap-cases.aliquots-graph.vertex.json.gz",
		"graph/graphbuild_depmap-cases.aliquots-graph.edge.json.gz",
		"graph/graphbuild_depmap-cases.depMapCases-graph.vertex.json.gz",
		"graph/graphbuild_depmap-cases.depMapCases-graph.edge.json.gz",
		"graph/graphbuild_depmap-cases.sampleObjects-graph.vertex.json.gz",
		"graph/graphbuild_depmap-cases.sampleObjects-graph.edge.json.gz"
	shell:
		"sifter run graph-build/graphbuild_depmap-cases.yaml"

rule graphbuild_depmap_expression:
	input:
		"output/depmap/depmap-expression.values.expression.json.gz",
		"schema",
		"graph-build/graphbuild_depmap-expression.yaml"
	output:
		"graph/graphbuild_depmap-expression.values-graph.edge.json.gz",
		"graph/graphbuild_depmap-expression.values-graph.vertex.json.gz"
	shell:
		"sifter run graph-build/graphbuild_depmap-expression.yaml"

rule graphbuild_depmap_mafs:
	input:
		"output/depmap/depmap-mafs.callsets.callset.json.gz",
		"output/depmap/depmap-mafs.variants.variants.json.gz",
		"schema",
		"schema",
		"graph-build/graphbuild_depmap-mafs.yaml"
	output:
		"graph/graphbuild_depmap-mafs.callsets-graph.vertex.json.gz",
		"graph/graphbuild_depmap-mafs.callsets-graph.edge.json.gz",
		"graph/graphbuild_depmap-mafs.variants-graph.vertex.json.gz",
		"graph/graphbuild_depmap-mafs.variants-graph.edge.json.gz"
	shell:
		"sifter run graph-build/graphbuild_depmap-mafs.yaml"

rule graphbuild_ensembl_gtf:
	input:
		"output/ensembl/ensembl_gtf.exons.exon.json.gz",
		"output/ensembl/ensembl_gtf.genes.gene.json.gz",
		"output/ensembl/ensembl_gtf.transcripts.transcript.json.gz",
		"schema",
		"schema",
		"schema",
		"graph-build/graphbuild_ensembl_gtf.yaml"
	output:
		"graph/graphbuild_ensembl_gtf.exons-graph.vertex.json.gz",
		"graph/graphbuild_ensembl_gtf.exons-graph.edge.json.gz",
		"graph/graphbuild_ensembl_gtf.genes-graph.vertex.json.gz",
		"graph/graphbuild_ensembl_gtf.genes-graph.edge.json.gz",
		"graph/graphbuild_ensembl_gtf.transcripts-graph.vertex.json.gz",
		"graph/graphbuild_ensembl_gtf.transcripts-graph.edge.json.gz"
	shell:
		"sifter run graph-build/graphbuild_ensembl_gtf.yaml"

rule graphbuild_g2p:
	input:
		"output/g2p/g2p.main.assocation.json.gz",
		"schema",
		"graph-build/graphbuild_g2p.yaml"
	output:
		"graph/graphbuild_g2p.main-graph.vertex.json.gz",
		"graph/graphbuild_g2p.main-graph.edge.json.gz"
	shell:
		"sifter run graph-build/graphbuild_g2p.yaml"

rule graphbuild_gdc_mafs:
	input:
		"output/gdc/open-maf/gdc-mafs.scan.variant.json.gz",
		"output/gdc/open-maf/gdc-mafs.somaticCallsets.callset.json.gz",
		"schema",
		"schema",
		"graph-build/graphbuild_gdc-mafs.yaml"
	output:
		"graph/graphbuild_gdc-mafs.somaticCallsets-graph.edge.json.gz",
		"graph/graphbuild_gdc-mafs.scan-graph.vertex.json.gz",
		"graph/graphbuild_gdc-mafs.scan-graph.edge.json.gz",
		"graph/graphbuild_gdc-mafs.somaticCallsets-graph.vertex.json.gz"
	shell:
		"sifter run graph-build/graphbuild_gdc-mafs.yaml"

rule graphbuild_gdc:
	input:
		"output/gdc/gdc.aliquotObject.aliquot.json.gz",
		"output/gdc/gdc.caseObject.case.json.gz",
		"output/gdc/gdc.projectObject.project.json.gz",
		"output/gdc/gdc.sampleObject.sample.json.gz",
		"schema",
		"schema",
		"schema",
		"schema",
		"graph-build/graphbuild_gdc.yaml"
	output:
		"graph/graphbuild_gdc.aliquotObject-graph.edge.json.gz",
		"graph/graphbuild_gdc.caseObject-graph.vertex.json.gz",
		"graph/graphbuild_gdc.caseObject-graph.edge.json.gz",
		"graph/graphbuild_gdc.projectObject-graph.vertex.json.gz",
		"graph/graphbuild_gdc.projectObject-graph.edge.json.gz",
		"graph/graphbuild_gdc.sampleObject-graph.vertex.json.gz",
		"graph/graphbuild_gdc.sampleObject-graph.edge.json.gz",
		"graph/graphbuild_gdc.aliquotObject-graph.vertex.json.gz"
	shell:
		"sifter run graph-build/graphbuild_gdc.yaml"

rule graphbuild_go:
	input:
		"output/go/go.transform.term.json.gz",
		"schema",
		"graph-build/graphbuild_go.yaml"
	output:
		"graph/graphbuild_go.transform-graph.vertex.json.gz",
		"graph/graphbuild_go.transform-graph.edge.json.gz"
	shell:
		"sifter run graph-build/graphbuild_go.yaml"

rule graphbuild_mondo:
	input:
		"output/mondo/mondo.extract.phenotype.json.gz",
		"schema",
		"graph-build/graphbuild_mondo.yaml"
	output:
		"graph/graphbuild_mondo.extract-graph.edge.json.gz",
		"graph/graphbuild_mondo.extract-graph.vertex.json.gz"
	shell:
		"sifter run graph-build/graphbuild_mondo.yaml"

rule graphbuild_msigdb:
	input:
		"output/msigdb/msigdb.transform.gene_set.json.gz",
		"schema",
		"graph-build/graphbuild_msigdb.yaml"
	output:
		"graph/graphbuild_msigdb.transform-graph.vertex.json.gz",
		"graph/graphbuild_msigdb.transform-graph.edge.json.gz"
	shell:
		"sifter run graph-build/graphbuild_msigdb.yaml"

rule graphbuild_pathway_commons:
	input:
		"output/pathway_commons/pathway_commons.interactionMap.interaction.json.gz",
		"schema",
		"graph-build/graphbuild_pathway_commons.yaml"
	output:
		"graph/graphbuild_pathway_commons.interactionMap-graph.vertex.json.gz",
		"graph/graphbuild_pathway_commons.interactionMap-graph.edge.json.gz"
	shell:
		"sifter run graph-build/graphbuild_pathway_commons.yaml"

rule graphbuild_pdb:
	input:
		"output/pdb/pdb.start.protein_structure.json.gz",
		"schema",
		"graph-build/graphbuild_pdb.yaml"
	output:
		"graph/graphbuild_pdb.start-graph.vertex.json.gz",
		"graph/graphbuild_pdb.start-graph.edge.json.gz"
	shell:
		"sifter run graph-build/graphbuild_pdb.yaml"

rule graphbuild_pharmacodb_profiles_gx:
	input:
		"output/pharmacodb/pharmacodb_profiles_gx.aliquot.aliquot.json.gz",
		"output/pharmacodb/pharmacodb_profiles_gx.drugResponse.drug_response.json.gz",
		"output/pharmacodb/pharmacodb_profiles_gx.sample.sample.json.gz",
		"schema",
		"schema",
		"schema",
		"graph-build/graphbuild_pharmacodb_profiles_gx.yaml"
	output:
		"graph/graphbuild_pharmacodb_profiles_gx.aliquot-graph.edge.json.gz",
		"graph/graphbuild_pharmacodb_profiles_gx.drugResponse-graph.vertex.json.gz",
		"graph/graphbuild_pharmacodb_profiles_gx.drugResponse-graph.edge.json.gz",
		"graph/graphbuild_pharmacodb_profiles_gx.sample-graph.vertex.json.gz",
		"graph/graphbuild_pharmacodb_profiles_gx.sample-graph.edge.json.gz",
		"graph/graphbuild_pharmacodb_profiles_gx.aliquot-graph.vertex.json.gz"
	shell:
		"sifter run graph-build/graphbuild_pharmacodb_profiles_gx.yaml"

rule graphbuild_prism_transform:
	input:
		"output/prism/prism_transform.secondary.drug_response.json.gz",
		"output/prism/prism_transform.primary.drug_response.json.gz",
		"schema",
		"schema",
		"graph-build/graphbuild_prism_transform.yaml"
	output:
		"graph/graphbuild_prism_transform.primary-graph.vertex.json.gz",
		"graph/graphbuild_prism_transform.primary-graph.edge.json.gz",
		"graph/graphbuild_prism_transform.secondary-graph.vertex.json.gz",
		"graph/graphbuild_prism_transform.secondary-graph.edge.json.gz"
	shell:
		"sifter run graph-build/graphbuild_prism_transform.yaml"

rule graphbuild_pubmed:
	input:
		"output/pubmed/pubmed.transform.publication.json.gz",
		"schema",
		"graph-build/graphbuild_pubmed.yaml"
	output:
		"graph/graphbuild_pubmed.transform-graph.vertex.json.gz",
		"graph/graphbuild_pubmed.transform-graph.edge.json.gz"
	shell:
		"sifter run graph-build/graphbuild_pubmed.yaml"

rule graphbuild_rnaseq:
	input:
		"output/gdc/rnaseq.rna.gene_rnaseq.json.gz",
		"schema",
		"graph-build/graphbuild_rnaseq.yaml"
	output:
		"graph/graphbuild_rnaseq.rna-graph.vertex.json.gz",
		"graph/graphbuild_rnaseq.rna-graph.edge.json.gz"
	shell:
		"sifter run graph-build/graphbuild_rnaseq.yaml"

rule graphbuild_uniprot_sprot:
	input:
		"output/uniprot/uniprot_sprot.start.protein.json.gz",
		"schema",
		"graph-build/graphbuild_uniprot_sprot.yaml"
	output:
		"graph/graphbuild_uniprot_sprot.start-graph.edge.json.gz",
		"graph/graphbuild_uniprot_sprot.start-graph.vertex.json.gz"
	shell:
		"sifter run graph-build/graphbuild_uniprot_sprot.yaml"

rule graphbuild_uniprot_trembl:
	input:
		"output/uniprot/uniprot_trembl.start.protein.json.gz",
		"schema",
		"graph-build/graphbuild_uniprot_trembl.yaml"
	output:
		"graph/graphbuild_uniprot_trembl.start-graph.vertex.json.gz",
		"graph/graphbuild_uniprot_trembl.start-graph.edge.json.gz"
	shell:
		"sifter run graph-build/graphbuild_uniprot_trembl.yaml"

