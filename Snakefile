


rule all:
	input:
		"source/pathway_commons/PathwayCommons12.panther.complex",
		"output/uniprot/uniprot_trembl.start.protein.json.gz",
		"output/bindingdb/bindingdb.assay.assay.json.gz",
		"output/chembl/chemblTransform.records.compound.json.gz",
		"source/pathway_commons/PathwayCommons12.corum.extSIF",
		"source/gdc/files.json",
		"output/gdc/rnaseq.rna.gene_rnaseq.json.gz",
		"output/mondo/mondo.extract.phenotype.json.gz",
		"source/pathway_commons/PathwayCommons12.panther.extSIF",
		"source/pathway_commons/PathwayCommons12.inoh.extSIF",
		"source/pathway_commons/PathwayCommons12.pid.extSIF",
		"output/pathway_commons/pathway_commons.complexBundle.complex.json.gz",
		"source/pathway_commons/PathwayCommons12.kegg.complex",
		"source/pathway_commons/PathwayCommons12.msigdb.complex",
		"source/pathway_commons/PathwayCommons12.msigdb.extSIF",
		"source/pathway_commons/PathwayCommons12.innatedb.extSIF",
		"output/prism/prism_transform.primary.drug_response.json.gz",
		"output/gdc/gdc.aliquotAlias.table.json.gz",
		"output/gtex/GTEX_Gene_Expression.gctProcess.gene_expression.json.gz",
		"source/pathway_commons/PathwayCommons12.netpath.complex",
		"source/pathway_commons/PathwayCommons12.biogrid.extSIF",
		"output/bindingdb/bindingdbTsv.start.protein_compound_association.json.gz",
		"output/gdsc/GDSC_Transform.aliquot.aliquot.json.gz",
		"output/bindingdb/bindingdb.start.raw.json.gz",
		"source/pharmacodb/tables/source_statistics.tsv.gz",
		"source/pharmacodb/tables/source_tissue_names.tsv.gz",
		"source/pharmacodb/tables/gene_drugs.tsv.gz",
		"output/pharmacodb/pharmacodb_profiles.cellAliquot.aliquot.json.gz",
		"output/bindingdb/bindingdbTsv.debug.debug.json.gz",
		"output/g2p/g2p.main.assocation.json.gz",
		"source/pathway_commons/PathwayCommons12.inoh.complex",
		"source/gtex/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt",
		"output/ensembl/ensembl_gtf.genes.gene.json.gz",
		"source/ensembl/Homo_sapiens.GRCh37.85.uniprot.tsv.gz",
		"source/pharmacodb/tables/source_cell_names.tsv.gz",
		"source/pathway_commons/PathwayCommons12.drugbank.complex",
		"output/pdb/pdb.start.protein_structure.json.gz",
		"source/pharmacodb/tables/dataset_cells.tsv.gz",
		"source/go/go.json",
		"source/pathway_commons/PathwayCommons12.netpath.extSIF",
		"output/uniprot/uniprot_sprot.start.protein.json.gz",
		"source/pathway_commons/PathwayCommons12.psp.complex",
		"output/chembl/chemblDrugMechanismExtract.build.protein_compound_association.json.gz",
		"output/ensembl/ensembl_gtf.exons.exon.json.gz",
		"output/gdc/gdc.aliquotObject.aliquot.json.gz",
		"source/pathway_commons/PathwayCommons12.humancyc.complex",
		"source/g2p/tables/g2p_source.bed",
		"output/pathway_commons/pathway_commons.interactionMap.interaction.json.gz",
		"output/msigdb/msigdb.transform.gene_set.json.gz",
		"source/pathway_commons/PathwayCommons12.reactome.complex",
		"source/pathway_commons/PathwayCommons12.hprd.complex",
		"source/pathway_commons/PathwayCommons12.bind.extSIF",
		"output/pharmacodb/pharmacodb_profiles.cellProject.project.json.gz",
		"source/pathway_commons/PathwayCommons12.psp.extSIF",
		"source/pharmacodb/tables/tissues.tsv.gz",
		"source/docm/variants.json",
		"output/gdc/gdc.caseObject.case.json.gz",
		"source/pathway_commons/PathwayCommons12.reactome.extSIF",
		"source/pharmgkb/relationships.zip",
		"output/bindingdb/bindingdb.monomer.monomer.json.gz",
		"source/pathway_commons/PathwayCommons12.mirtarbase.extSIF",
		"source/pathway_commons/PathwayCommons12.pathbank.extSIF",
		"source/pathway_commons/PathwayCommons12.bind.complex",
		"source/pathway_commons/PathwayCommons12.pathbank.complex",
		"output/depmap/depmap-expression.values.expression.json.gz",
		"output/gdc/gdc.sampleObject.sample.json.gz",
		"output/go/go.transform.term.json.gz",
		"source/pathway_commons/PathwayCommons12.ctd.extSIF",
		"source/ucscGenome/cytoBandIdeo.txt.gz",
		"output/gdsc/GDSC_rnaseq_Transform.start.geneExpression.json.gz",
		"source/pathway_commons/PathwayCommons12.reconx.extSIF",
		"source/pathway_commons/PathwayCommons12.hprd.extSIF",
		"source/pathway_commons/PathwayCommons12.innatedb.complex",
		"output/ensembl/ensembl_gtf.transcripts.transcript.json.gz",
		"output/gdsc/GDSC_rnaseq_Transform.aliquot.aliquot.json.gz",
		"source/pathway_commons/PathwayCommons12.ctd.complex",
		"output/depmap/depmap-cases.aliquots.aliquot.json.gz",
		"output/gdsc/GDSC_VCF_Transform.variants.somatic_variant.json.gz",
		"output/gtex/GTEX_Transcript_Expression.update.transcript_expression.json.gz",
		"output/depmap/depmap-mafs.callsets.callset.json.gz",
		"source/pharmacodb/tables/cell_tissues.tsv.gz",
		"output/pubmed/pubmed.transform.publication.json.gz",
		"source/pharmacodb/tables/drug_annots.tsv.gz",
		"source/pharmacodb/tables/cellosaurus.tsv.gz",
		"output/gdsc/GDSC_Transform.transform.geneExpression.json.gz",
		"source/ncit/ncit.obo",
		"source/pathway_commons/PathwayCommons12.reconx.complex",
		"source/pathway_commons/PathwayCommons12.dip.extSIF",
		"source/pathway_commons/PathwayCommons12.humancyc.extSIF",
		"output/pharmacodb/pharmacodb_profiles.drObject.drug_response.json.gz",
		"source/pharmgkb/clinicalAnnotations.zip",
		"output/gdsc/GDSC_VCF_Transform.callset.somatic_callset.json.gz",
		"source/pathway_commons/PathwayCommons12.biogrid.complex",
		"source/pharmacodb/tables/source_drug_names.tsv.gz",
		"output/pharmacodb/pharmacodb_profiles.cellDistinct.checkpoint.json.gz",
		"output/cellosarus/cellosarus.transform.cases.json.gz",
		"output/gdc/open-maf/gdc-mafs.somaticCallsets.callset.json.gz",
		"source/pharmacodb/tables/sources.tsv.gz",
		"output/gdc/open-maf/gdc-mafs.scan.variant.json.gz",
		"source/pathway_commons/PathwayCommons12.kegg.extSIF",
		"source/pathway_commons/PathwayCommons12.mirtarbase.complex",
		"source/pathway_commons/PathwayCommons12.pid.complex",
		"source/pathway_commons/PathwayCommons12.corum.complex",
		"output/prism/prism_transform.secondary.drug_response.json.gz",
		"output/clinvar/clinvar.transform.raw.json.gz",
		"output/depmap/depmap-mafs.variants.variants.json.gz",
		"output/gdsc/GDSC_VCF_Transform.namefix.source.json.gz",
		"source/pathway_commons/PathwayCommons12.dip.complex",
		"output/allele/annotatedAllele.transform.allele.json.gz",
		"output/depmap/depmap-cases.samples.sample.json.gz",
		"output/pharmacodb/pharmacodb_profiles.cellSample.sample.json.gz",
		"source/gtex/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt",
		"source/pathway_commons/PathwayCommons12.drugbank.extSIF",
		"output/pharmacodb/pharmacodb_profiles.cellCase.case.json.gz",
		"output/bindingdb/bindingdb.complex.complex.json.gz",
		"source/ensembl/Homo_sapiens.GRCh38.108.uniprot.tsv.gz",
		"source/ensembl/homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20221007.gff.gz",
		"output/go/go_gaf.dump.gaf.json.gz"

rule bindingdb:
	input:
		"source/bindingdb/BDB_my-202301.dmp",
		"transform/bindingdb/transform.yaml"
	output:
		"output/bindingdb/bindingdb.assay.assay.json.gz",
		"output/bindingdb/bindingdb.complex.complex.json.gz",
		"output/bindingdb/bindingdb.monomer.monomer.json.gz",
		"output/bindingdb/bindingdb.start.raw.json.gz"
	shell:
		"sifter run transform/bindingdb/transform.yaml"

rule bindingdbTsv:
	input:
		"source/bindingdb/BindingDB_All.tsv",
		"schema",
		"transform/bindingdb/transform_tsv.yaml"
	output:
		"output/bindingdb/bindingdbTsv.debug.debug.json.gz",
		"output/bindingdb/bindingdbTsv.start.protein_compound_association.json.gz"
	shell:
		"sifter run transform/bindingdb/transform_tsv.yaml"

rule synonyms:
	input:
		"source/cellosarus/cellosaurus.obo",
		"transform/cellosarus/alias_table.yaml"
	output:
		"source/cellosarus/table/synonyms.caseTable.ach2cellosarus.json.gz"
	shell:
		"sifter run transform/cellosarus/alias_table.yaml"

rule cellosarus:
	input:
		"source/cellosarus/cellosaurus.obo",
		"source/mondo/tables/mondo.ncit_extract.mapping.json.gz",
		"schema",
		"transform/cellosarus/transform.yaml"
	output:
		"output/cellosarus/cellosarus.transform.cases.json.gz"
	shell:
		"sifter run transform/cellosarus/transform.yaml"

rule chemblDownload_curl:
	output:
		"source/chembl/chembl_31_sqlite.tar.gz"
	shell:
		"cd transform/chembl && curl -o ../../source/chembl/chembl_31_sqlite.tar.gz https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/chembl_31_sqlite.tar.gz"

rule chemblDownload_tar:
	input:
		"source/chembl/chembl_31_sqlite.tar.gz"
	output:
		"source/chembl/chembl_31/chembl_31_sqlite/chembl_31.db"
	shell:
		"cd transform/chembl && tar xvzf ../../source/chembl/chembl_31_sqlite.tar.gz -C ../../source/chembl/"

rule chemblDrugMechanismExtract:
	input:
		"source/chembl/chembl_31/chembl_31_sqlite/chembl_31.db",
		"schema",
		"transform/chembl/drug_mechanism.yaml"
	output:
		"output/chembl/chemblDrugMechanismExtract.build.protein_compound_association.json.gz"
	shell:
		"sifter run transform/chembl/drug_mechanism.yaml"

rule synonyms_1:
	input:
		"source/chembl/chembl_31/chembl_31_sqlite/chembl_31.db",
		"transform/chembl/synonyms.yaml"
	output:
		"source/chembl/tables/synonyms.buildTable.synonyms.json.gz"
	shell:
		"sifter run transform/chembl/synonyms.yaml"

rule chemblTransform:
	input:
		"source/chembl/chembl_31/chembl_31_sqlite/chembl_31.db",
		"source/chembl/tables/synonyms.buildTable.synonyms.json.gz",
		"schema",
		"transform/chembl/transform.yaml"
	output:
		"output/chembl/chemblTransform.records.compound.json.gz"
	shell:
		"sifter run transform/chembl/transform.yaml"

rule clinvar:
	input:
		"source/clinvar/ClinVarFullRelease_2023-01.xml.gz",
		"transform/clinvar/transform.yaml"
	output:
		"output/clinvar/clinvar.transform.raw.json.gz"
	shell:
		"sifter run transform/clinvar/transform.yaml"

rule depmap_cases:
	input:
		"source/depmap/Model.csv",
		"schema",
		"source/cellosarus/table/synonyms.caseTable.ach2cellosarus.json.gz",
		"schema",
		"transform/depmap/cases.yaml"
	output:
		"output/depmap/depmap-cases.samples.sample.json.gz",
		"output/depmap/depmap-cases.aliquots.aliquot.json.gz"
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
		"schema",
		"transform/depmap/mutations.yaml"
	output:
		"output/depmap/depmap-mafs.allele.allele.json.gz",
		"output/depmap/depmap-mafs.callsets.callset.json.gz",
		"output/depmap/depmap-mafs.variants.variants.json.gz"
	shell:
		"sifter run transform/depmap/mutations.yaml"

rule download_downloadGO:
	output:
		"source/docm/variants.json"
	shell:
		"cd transform/docm && curl -o ../../source/docm/variants.json http://www.docm.info/api/v1/variants.json"

rule ensembleDownload_ensemblDownload:
	output:
		"source/ensembl/Homo_sapiens.GRCh38.108.chr_patch_hapl_scaff.gtf.gz"
	shell:
		"cd transform/ensembl && curl -o ../../source/ensembl/Homo_sapiens.GRCh38.108.chr_patch_hapl_scaff.gtf.gz https://ftp.ensembl.org/pub/release-108/gtf/homo_sapiens/Homo_sapiens.GRCh38.108.chr_patch_hapl_scaff.gtf.gz"

rule ensembleDownload_ensemblProteinDownload:
	output:
		"source/ensembl/Homo_sapiens.GRCh38.108.uniprot.tsv.gz"
	shell:
		"cd transform/ensembl && curl -o ../../source/ensembl/Homo_sapiens.GRCh38.108.uniprot.tsv.gz https://ftp.ensembl.org/pub/release-108/tsv/homo_sapiens/Homo_sapiens.GRCh38.108.uniprot.tsv.gz"

rule ensembleDownload_regulatory:
	output:
		"source/ensembl/homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20221007.gff.gz"
	shell:
		"cd transform/ensembl && curl -o ../../source/ensembl/homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20221007.gff.gz  https://ftp.ensembl.org/pub/current_regulation/homo_sapiens/homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20221007.gff.gz"

rule ensembl_gtf:
	input:
		"source/ensembl/Homo_sapiens.GRCh38.108.chr_patch_hapl_scaff.gtf.gz",
		"schema",
		"schema",
		"schema",
		"transform/ensembl/ensembl_transform.yaml"
	output:
		"output/ensembl/ensembl_gtf.exons.exon.json.gz",
		"output/ensembl/ensembl_gtf.genes.gene.json.gz",
		"output/ensembl/ensembl_gtf.transcripts.transcript.json.gz"
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

rule g2p:
	input:
		"source/g2p",
		"transform/g2p/bed_file.yaml"
	output:
		"source/g2p/tables/g2p_source.bed"
	shell:
		"sifter run transform/g2p/bed_file.yaml"

rule g2p_1:
	input:
		"source/g2p",
		"source/hugo/hugo.tsv",
		"source/g2p/tables/hglft_genome_2749d_8437f0.bed",
		"schema",
		"schema",
		"transform/g2p/transform.yaml"
	output:
		"output/g2p/g2p.main.assocation.json.gz",
		"output/g2p/g2p.alleles.allele.json.gz"
	shell:
		"sifter run transform/g2p/transform.yaml"

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
		"schema",
		"transform/gdc/maf-files.yaml"
	output:
		"output/gdc/open-maf/gdc-mafs.somaticCallsets.callset.json.gz",
		"output/gdc/open-maf/gdc-mafs.allele.allele.json.gz",
		"output/gdc/open-maf/gdc-mafs.scan.variant.json.gz"
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

rule cosmic2ach:
	input:
		"source/gdsc/model_list_20191104.csv",
		"transform/gdsc/cosmic_to_ach.yaml"
	output:
		"source/gdsc/tables/cosmic2ach.translate.link.json.gz"
	shell:
		"sifter run transform/gdsc/cosmic_to_ach.yaml"

rule downloadGDSC_geneInfo:
	output:
		"source/gdsc/gene_identifiers_20191101.csv"
	shell:
		"cd transform/gdsc && curl https://cog.sanger.ac.uk/cmp/download/gene_identifiers_20191101.csv -o ../../source/gdsc/gene_identifiers_20191101.csv"

rule downloadGDSC_rnaSepUnzip:
	input:
		"source/gdsc/rnaseq_sanger_20210316.zip"
	output:
		"source/gdsc/rnaseq_sanger_20210316.csv"
	shell:
		"cd transform/gdsc && unzip -d ../../source/gdsc ../../source/gdsc/rnaseq_sanger_20210316.zip"

rule downloadGDSC_rnaSeq:
	output:
		"source/gdsc/rnaseq_sanger_20210316.zip"
	shell:
		"cd transform/gdsc && curl https://cog.sanger.ac.uk/cmp/download/rnaseq_sanger_20210316.zip -o ../../source/gdsc/rnaseq_sanger_20210316.zip"

rule downloadGDSC_sampleInfo:
	output:
		"source/gdsc/model_list_20230110.csv"
	shell:
		"cd transform/gdsc && curl https://cog.sanger.ac.uk/cmp/download/model_list_20230110.csv -o ../../source/gdsc/model_list_20230110.csv"

rule GDSC_Transform:
	input:
		"source/gdsc/Cell_line_RMA_proc_basalExp.txt",
		"source/gdsc/tables/cosmic2ach.translate.link.json.gz",
		"schema",
		"source/hugo/hugo.tsv",
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
		"output/gdsc/GDSC_VCF_Transform.variants.somatic_variant.json.gz",
		"output/gdsc/GDSC_VCF_Transform.allele.allele.json.gz",
		"output/gdsc/GDSC_VCF_Transform.callset.somatic_callset.json.gz",
		"output/gdsc/GDSC_VCF_Transform.namefix.source.json.gz"
	shell:
		"sifter run transform/gdsc/vcf_transform.yaml"

rule download_downloadGAF:
	output:
		"source/go/goa_human.gaf.gz"
	shell:
		"cd transform/go && curl -o ../../source/go/goa_human.gaf.gz http://release.geneontology.org/2022-09-19/annotations/goa_human.gaf.gz"

rule download_downloadGO_1:
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

rule go:
	input:
		"source/go/go.obo",
		"schema",
		"transform/go/go_obo.yaml"
	output:
		"output/go/go.transform.term.json.gz"
	shell:
		"sifter run transform/go/go_obo.yaml"

rule download_gtex_downloadSubjectPhenotypes:
	output:
		"source/gtex/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt"
	shell:
		"cd transform/gtex && curl -o ../../source/gtex/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt"

rule download_gtex_downloadTranscriptTPM:
	output:
		"source/gtex/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct.gz"
	shell:
		"cd transform/gtex && curl https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct.gz -o ../../source/gtex/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct.gz"

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

rule GTEX_Gene_Expression:
	input:
		"source/gtex/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz",
		"schema",
		"transform/gtex/gene_transform.yaml"
	output:
		"output/gtex/GTEX_Gene_Expression.gctProcess.gene_expression.json.gz"
	shell:
		"sifter run transform/gtex/gene_transform.yaml"

rule GTEX_Transcript_Expression:
	input:
		"source/gtex/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct.gz",
		"schema",
		"transform/gtex/transcript_transform.yaml"
	output:
		"output/gtex/GTEX_Transcript_Expression.update.transcript_expression.json.gz"
	shell:
		"sifter run transform/gtex/transcript_transform.yaml"

rule mondo:
	input:
		"source/mondo/mondo.json",
		"schema",
		"transform/mondo/mondo.yaml"
	output:
		"output/mondo/mondo.extract.phenotype.json.gz"
	shell:
		"sifter run transform/mondo/mondo.yaml"

rule mondo_1:
	input:
		"source/mondo/mondo.json",
		"transform/mondo/ncit2mondo_table.yaml"
	output:
		"source/mondo/tables/mondo.ncit_extract.mapping.json.gz"
	shell:
		"sifter run transform/mondo/ncit2mondo_table.yaml"

rule downloads_msigdb:
	output:
		"source/msigdb/msigdb_v7.5.1.xml"
	shell:
		"cd transform/msigdb && curl -o ../../source/msigdb/msigdb_v7.5.1.xml https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.5.1/msigdb_v7.5.1.xml"

rule msigdb:
	input:
		"source/msigdb/msigdb_v7.5.1.xml",
		"source/ensembl/gene2ensembl.translate.link.json.gz",
		"schema",
		"transform/msigdb/transform.yaml"
	output:
		"output/msigdb/msigdb.transform.gene_set.json.gz"
	shell:
		"sifter run transform/msigdb/transform.yaml"

rule download_downloadGO_2:
	output:
		"source/ncit/ncit.obo"
	shell:
		"cd transform/ncit && curl -L -o ../../source/ncit/ncit.obo https://github.com/NCI-Thesaurus/thesaurus-obo-edition/releases/download/v2022-08-19/ncit.obo"

rule pathwayCommonsDownload_pathwayCommonsDownload_pc_prep_download_psp:
	output:
		"source/pathway_commons/PathwayCommons12.psp.BIOPAX.owl.gz"
	shell:
		"cd transform/pathway_commons && curl -o ../../source/pathway_commons/PathwayCommons12.psp.BIOPAX.owl.gz https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.psp.BIOPAX.owl.gz"

rule pathwayCommonsDownload_pathwayCommonsDownload_pc_prep_download_msigdb:
	output:
		"source/pathway_commons/PathwayCommons12.msigdb.BIOPAX.owl.gz"
	shell:
		"cd transform/pathway_commons && curl -o ../../source/pathway_commons/PathwayCommons12.msigdb.BIOPAX.owl.gz https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.msigdb.BIOPAX.owl.gz"

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
		"cd transform/pathway_commons && java -Xmx15g -jar ../../tools/extsif/target/pc-extract-1.0-jar-with-dependencies.jar ../../source/pathway_commons/PathwayCommons12.kegg.BIOPAX.owl.gz ../../source/pathway_commons/PathwayCommons12.kegg"

rule pathwayCommonsDownload_pathwayCommonsDownload_pc_prep_download_hprd:
	output:
		"source/pathway_commons/PathwayCommons12.hprd.BIOPAX.owl.gz"
	shell:
		"cd transform/pathway_commons && curl -o ../../source/pathway_commons/PathwayCommons12.hprd.BIOPAX.owl.gz https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.hprd.BIOPAX.owl.gz"

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
		"cd transform/pathway_commons && java -Xmx15g -jar ../../tools/extsif/target/pc-extract-1.0-jar-with-dependencies.jar ../../source/pathway_commons/PathwayCommons12.ctd.BIOPAX.owl.gz ../../source/pathway_commons/PathwayCommons12.ctd"

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
		"cd transform/pathway_commons && java -Xmx15g -jar ../../tools/extsif/target/pc-extract-1.0-jar-with-dependencies.jar ../../source/pathway_commons/PathwayCommons12.humancyc.BIOPAX.owl.gz ../../source/pathway_commons/PathwayCommons12.humancyc"

rule pathwayCommonsDownload_pathwayCommonsDownload_pc_prep_download_drugbank:
	output:
		"source/pathway_commons/PathwayCommons12.drugbank.BIOPAX.owl.gz"
	shell:
		"cd transform/pathway_commons && curl -o ../../source/pathway_commons/PathwayCommons12.drugbank.BIOPAX.owl.gz https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.drugbank.BIOPAX.owl.gz"

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
		"cd transform/pathway_commons && java -Xmx15g -jar ../../tools/extsif/target/pc-extract-1.0-jar-with-dependencies.jar ../../source/pathway_commons/PathwayCommons12.mirtarbase.BIOPAX.owl.gz ../../source/pathway_commons/PathwayCommons12.mirtarbase"

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
		"cd transform/pathway_commons && java -Xmx15g -jar ../../tools/extsif/target/pc-extract-1.0-jar-with-dependencies.jar ../../source/pathway_commons/PathwayCommons12.msigdb.BIOPAX.owl.gz ../../source/pathway_commons/PathwayCommons12.msigdb"

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
		"cd transform/pathway_commons && java -Xmx15g -jar ../../tools/extsif/target/pc-extract-1.0-jar-with-dependencies.jar ../../source/pathway_commons/PathwayCommons12.reconx.BIOPAX.owl.gz ../../source/pathway_commons/PathwayCommons12.reconx"

rule pathwayCommonsDownload_pathwayCommonsDownload_pc_prep_download_pathbank:
	output:
		"source/pathway_commons/PathwayCommons12.pathbank.BIOPAX.owl.gz"
	shell:
		"cd transform/pathway_commons && curl -o ../../source/pathway_commons/PathwayCommons12.pathbank.BIOPAX.owl.gz https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.pathbank.BIOPAX.owl.gz"

rule pathwayCommonsDownload_pathwayCommonsDownload_pc_prep_download_bind:
	output:
		"source/pathway_commons/PathwayCommons12.bind.BIOPAX.owl.gz"
	shell:
		"cd transform/pathway_commons && curl -o ../../source/pathway_commons/PathwayCommons12.bind.BIOPAX.owl.gz https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.bind.BIOPAX.owl.gz"

rule pathwayCommonsDownload_pathwayCommonsDownload_pc_prep_download_dip:
	output:
		"source/pathway_commons/PathwayCommons12.dip.BIOPAX.owl.gz"
	shell:
		"cd transform/pathway_commons && curl -o ../../source/pathway_commons/PathwayCommons12.dip.BIOPAX.owl.gz https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.dip.BIOPAX.owl.gz"

rule pathwayCommonsDownload_pathwayCommonsDownload_pc_prep_download_humancyc:
	output:
		"source/pathway_commons/PathwayCommons12.humancyc.BIOPAX.owl.gz"
	shell:
		"cd transform/pathway_commons && curl -o ../../source/pathway_commons/PathwayCommons12.humancyc.BIOPAX.owl.gz https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.humancyc.BIOPAX.owl.gz"

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
		"cd transform/pathway_commons && java -Xmx15g -jar ../../tools/extsif/target/pc-extract-1.0-jar-with-dependencies.jar ../../source/pathway_commons/PathwayCommons12.dip.BIOPAX.owl.gz ../../source/pathway_commons/PathwayCommons12.dip"

rule pathwayCommonsDownload_pathwayCommonsDownload_pc_prep_download_innatedb:
	output:
		"source/pathway_commons/PathwayCommons12.innatedb.BIOPAX.owl.gz"
	shell:
		"cd transform/pathway_commons && curl -o ../../source/pathway_commons/PathwayCommons12.innatedb.BIOPAX.owl.gz https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.innatedb.BIOPAX.owl.gz"

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
		"cd transform/pathway_commons && java -Xmx15g -jar ../../tools/extsif/target/pc-extract-1.0-jar-with-dependencies.jar ../../source/pathway_commons/PathwayCommons12.biogrid.BIOPAX.owl.gz ../../source/pathway_commons/PathwayCommons12.biogrid"

rule pathwayCommonsDownload_pathwayCommonsDownload_pc_prep_download_panther:
	output:
		"source/pathway_commons/PathwayCommons12.panther.BIOPAX.owl.gz"
	shell:
		"cd transform/pathway_commons && curl -o ../../source/pathway_commons/PathwayCommons12.panther.BIOPAX.owl.gz https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.panther.BIOPAX.owl.gz"

rule pathwayCommonsDownload_pathwayCommonsDownload_pc_prep_download_ctd:
	output:
		"source/pathway_commons/PathwayCommons12.ctd.BIOPAX.owl.gz"
	shell:
		"cd transform/pathway_commons && curl -o ../../source/pathway_commons/PathwayCommons12.ctd.BIOPAX.owl.gz https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.ctd.BIOPAX.owl.gz"

rule pathwayCommonsDownload_pathwayCommonsDownload_pc_prep_download_kegg:
	output:
		"source/pathway_commons/PathwayCommons12.kegg.BIOPAX.owl.gz"
	shell:
		"cd transform/pathway_commons && curl -o ../../source/pathway_commons/PathwayCommons12.kegg.BIOPAX.owl.gz https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.kegg.BIOPAX.owl.gz"

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
		"cd transform/pathway_commons && java -Xmx15g -jar ../../tools/extsif/target/pc-extract-1.0-jar-with-dependencies.jar ../../source/pathway_commons/PathwayCommons12.corum.BIOPAX.owl.gz ../../source/pathway_commons/PathwayCommons12.corum"

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
		"cd transform/pathway_commons && java -Xmx15g -jar ../../tools/extsif/target/pc-extract-1.0-jar-with-dependencies.jar ../../source/pathway_commons/PathwayCommons12.psp.BIOPAX.owl.gz ../../source/pathway_commons/PathwayCommons12.psp"

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
		"cd transform/pathway_commons && java -Xmx15g -jar ../../tools/extsif/target/pc-extract-1.0-jar-with-dependencies.jar ../../source/pathway_commons/PathwayCommons12.hprd.BIOPAX.owl.gz ../../source/pathway_commons/PathwayCommons12.hprd"

rule pathwayCommonsDownload_pathwayCommonsDownload_pc_prep_download_corum:
	output:
		"source/pathway_commons/PathwayCommons12.corum.BIOPAX.owl.gz"
	shell:
		"cd transform/pathway_commons && curl -o ../../source/pathway_commons/PathwayCommons12.corum.BIOPAX.owl.gz https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.corum.BIOPAX.owl.gz"

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
		"cd transform/pathway_commons && java -Xmx15g -jar ../../tools/extsif/target/pc-extract-1.0-jar-with-dependencies.jar ../../source/pathway_commons/PathwayCommons12.drugbank.BIOPAX.owl.gz ../../source/pathway_commons/PathwayCommons12.drugbank"

rule pathwayCommonsDownload_pathwayCommonsDownload_pc_prep_download_biogrid:
	output:
		"source/pathway_commons/PathwayCommons12.biogrid.BIOPAX.owl.gz"
	shell:
		"cd transform/pathway_commons && curl -o ../../source/pathway_commons/PathwayCommons12.biogrid.BIOPAX.owl.gz https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.biogrid.BIOPAX.owl.gz"

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
		"cd transform/pathway_commons && java -Xmx15g -jar ../../tools/extsif/target/pc-extract-1.0-jar-with-dependencies.jar ../../source/pathway_commons/PathwayCommons12.inoh.BIOPAX.owl.gz ../../source/pathway_commons/PathwayCommons12.inoh"

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
		"cd transform/pathway_commons && java -Xmx15g -jar ../../tools/extsif/target/pc-extract-1.0-jar-with-dependencies.jar ../../source/pathway_commons/PathwayCommons12.netpath.BIOPAX.owl.gz ../../source/pathway_commons/PathwayCommons12.netpath"

rule pathwayCommonsDownload_pathwayCommonsDownload_pc_prep_download_mirtarbase:
	output:
		"source/pathway_commons/PathwayCommons12.mirtarbase.BIOPAX.owl.gz"
	shell:
		"cd transform/pathway_commons && curl -o ../../source/pathway_commons/PathwayCommons12.mirtarbase.BIOPAX.owl.gz https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.mirtarbase.BIOPAX.owl.gz"

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
		"cd transform/pathway_commons && java -Xmx15g -jar ../../tools/extsif/target/pc-extract-1.0-jar-with-dependencies.jar ../../source/pathway_commons/PathwayCommons12.pid.BIOPAX.owl.gz ../../source/pathway_commons/PathwayCommons12.pid"

rule pathwayCommonsDownload_pathwayCommonsDownload_pc_prep_download_reconx:
	output:
		"source/pathway_commons/PathwayCommons12.reconx.BIOPAX.owl.gz"
	shell:
		"cd transform/pathway_commons && curl -o ../../source/pathway_commons/PathwayCommons12.reconx.BIOPAX.owl.gz https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.reconx.BIOPAX.owl.gz"

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
		"cd transform/pathway_commons && java -Xmx15g -jar ../../tools/extsif/target/pc-extract-1.0-jar-with-dependencies.jar ../../source/pathway_commons/PathwayCommons12.bind.BIOPAX.owl.gz ../../source/pathway_commons/PathwayCommons12.bind"

rule pathwayCommonsDownload_pathwayCommonsDownload_pc_prep_download_inoh:
	output:
		"source/pathway_commons/PathwayCommons12.inoh.BIOPAX.owl.gz"
	shell:
		"cd transform/pathway_commons && curl -o ../../source/pathway_commons/PathwayCommons12.inoh.BIOPAX.owl.gz https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.inoh.BIOPAX.owl.gz"

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
		"cd transform/pathway_commons && java -Xmx15g -jar ../../tools/extsif/target/pc-extract-1.0-jar-with-dependencies.jar ../../source/pathway_commons/PathwayCommons12.innatedb.BIOPAX.owl.gz ../../source/pathway_commons/PathwayCommons12.innatedb"

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
		"cd transform/pathway_commons && java -Xmx15g -jar ../../tools/extsif/target/pc-extract-1.0-jar-with-dependencies.jar ../../source/pathway_commons/PathwayCommons12.reactome.BIOPAX.owl.gz ../../source/pathway_commons/PathwayCommons12.reactome"

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
		"cd transform/pathway_commons && java -Xmx15g -jar ../../tools/extsif/target/pc-extract-1.0-jar-with-dependencies.jar ../../source/pathway_commons/PathwayCommons12.panther.BIOPAX.owl.gz ../../source/pathway_commons/PathwayCommons12.panther"

rule pathwayCommonsDownload_pathwayCommonsDownload_pc_prep_download_reactome:
	output:
		"source/pathway_commons/PathwayCommons12.reactome.BIOPAX.owl.gz"
	shell:
		"cd transform/pathway_commons && curl -o ../../source/pathway_commons/PathwayCommons12.reactome.BIOPAX.owl.gz https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.reactome.BIOPAX.owl.gz"

rule pathwayCommonsDownload_pathwayCommonsDownload_pc_prep_download_netpath:
	output:
		"source/pathway_commons/PathwayCommons12.netpath.BIOPAX.owl.gz"
	shell:
		"cd transform/pathway_commons && curl -o ../../source/pathway_commons/PathwayCommons12.netpath.BIOPAX.owl.gz https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.netpath.BIOPAX.owl.gz"

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
		"cd transform/pathway_commons && java -Xmx15g -jar ../../tools/extsif/target/pc-extract-1.0-jar-with-dependencies.jar ../../source/pathway_commons/PathwayCommons12.pathbank.BIOPAX.owl.gz ../../source/pathway_commons/PathwayCommons12.pathbank"

rule pathwayCommonsDownload_pathwayCommonsDownload_pc_prep_download_pid:
	output:
		"source/pathway_commons/PathwayCommons12.pid.BIOPAX.owl.gz"
	shell:
		"cd transform/pathway_commons && curl -o ../../source/pathway_commons/PathwayCommons12.pid.BIOPAX.owl.gz https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.pid.BIOPAX.owl.gz"

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

rule pdb:
	input:
		"source/pdb/entries.idx",
		"schema",
		"transform/pdb/transform.yaml"
	output:
		"output/pdb/pdb.start.protein_structure.json.gz"
	shell:
		"sifter run transform/pdb/transform.yaml"

rule pharmacodb:
	input:
		"source/pharmacodb/pharmacodb-1.1.1.sql",
		"transform/pharmacodb/database.yaml"
	output:
		"source/pharmacodb/tables/drug_annots.tsv.gz",
		"source/pharmacodb/tables/dose_responses.tsv.gz",
		"source/pharmacodb/tables/profiles.tsv.gz",
		"source/pharmacodb/tables/source_statistics.tsv.gz",
		"source/pharmacodb/tables/source_tissue_names.tsv.gz",
		"source/pharmacodb/tables/gene_drugs.tsv.gz",
		"source/pharmacodb/tables/dataset_cells.tsv.gz",
		"source/pharmacodb/tables/drugs.tsv.gz",
		"source/pharmacodb/tables/source_cell_names.tsv.gz",
		"source/pharmacodb/tables/cell_tissues.tsv.gz",
		"source/pharmacodb/tables/cells.tsv.gz",
		"source/pharmacodb/tables/source_drug_names.tsv.gz",
		"source/pharmacodb/tables/cellosaurus.tsv.gz",
		"source/pharmacodb/tables/datasets.tsv.gz",
		"source/pharmacodb/tables/experiments.tsv.gz",
		"source/pharmacodb/tables/tissues.tsv.gz",
		"source/pharmacodb/tables/sources.tsv.gz"
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
		"source/pharmacodb/tables/dr_reduce.curveReduce.dose_response_curve.json.gz",
		"source/pharmacodb/tables/experiments.tsv.gz",
		"source/pharmacodb/tables/drugs.tsv.gz",
		"source/pharmacodb/tables/datasets.tsv.gz",
		"source/pharmacodb/tables/cells.tsv.gz",
		"schema",
		"schema",
		"schema",
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

rule downloadPharmgkb_downloadAnnotations:
	output:
		"source/pharmgkb/clinicalAnnotations.zip"
	shell:
		"cd transform/pharmgkb && curl -o ../../source/pharmgkb/clinicalAnnotations.zip https://api.pharmgkb.org/v1/download/file/data/clinicalAnnotations.zip"

rule downloadPharmgkb_downloadPharmgkb:
	output:
		"source/pharmgkb/relationships.zip"
	shell:
		"cd transform/pharmgkb && curl -o ../../source/pharmgkb/relationships.zip https://api.pharmgkb.org/v1/download/file/data/relationships.zip"

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
	shell:
		"sifter run transform/pubmed/transform.yaml"

rule downloadUCSC_downloadPharmgkb:
	output:
		"source/ucscGenome/cytoBandIdeo.txt.gz"
	shell:
		"cd transform/ucscGenome && curl -o ../../source/ucscGenome/cytoBandIdeo.txt.gz https://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/cytoBandIdeo.txt.gz"

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
		"normalize/allele/annotate_transform.yaml"
	output:
		"output/allele/annotatedAllele.transform.allele.json.gz"
	shell:
		"sifter run normalize/allele/annotate_transform.yaml"

rule alleleAnnotate_alleleAnnotate:
	input:
		"output-normalize/allele.vcf"
	output:
		"output-normalize/allele.annotated.vcf"
	shell:
		"cd normalize/allele && java -jar ../../util/snpEff/snpEff.jar ann -dataDir `pwd`/../../source/allele/data -nodownload  GRCh38.86 ../../output-normalize/allele.vcf > ../../output-normalize/allele.annotated.vcf"

rule annotated_VCF_transform:
	input:
		"output-normalize/allele.merge.json.gz",
		"normalize/allele/vcf_transform.yaml"
	output:
		"output-normalize/allele.vcf"
	shell:
		"sifter run normalize/allele/vcf_transform.yaml"



