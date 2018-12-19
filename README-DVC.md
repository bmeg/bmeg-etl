
DVC
-----

See:

* [Configure DVC](https://dvc.org/doc/get-started/configure)
* AWS CLI with CEPH

```
Once  you've got the awscli package installed and the "aws" utility in your PATH, you'll want to do two things to get S3 operations to work.

First, create a section in your ~/.aws/credentials file to hold your user keys (which will be sent to you individually in a separate message). In the example below, the profile will be named "ceph."

[ceph]
aws_access_key_id = 12345678901234567890
aws_secret_access_key = 1234567890123456789012345678901234567890

Second, and this is optional but I suggest it, create a shell alias. In my example, I call it "cephos" (Ceph Object Store). Make sure the profile argument matches the name in your AWS credentials file.

alias cephos="aws --profile=ceph --endpoint=https://ceph.acc.ohsu.edu"

```


Note: All dvc files are maintained at the project root.  This is due to:

 * each transformer is it's own command, executed with CWD at the project root
 * `dvc run` captures state in a separate \*.dvc file
 * therefore, there are quite a few files at the root level of the project




Setup
-----

```


# install and configure aws, see above
# test
$ cephos s3 ls --recursive s3://bmeg/dvc  | wc -l
654

# dvc already installed and initialized
# config should be ...
$ cat .dvc/config
# aws default credentials should point at ceph
['remote "ceph"']
url = s3://bmeg/dvc
endpointurl = https://ceph.acc.ohsu.edu
[core]
remote = ceph
['remote "source"']
url = s3://bmeg/source
endpointurl = https://ceph.acc.ohsu.edu
```

Example
----------

```
# retrieve data from foreign source
dvc run --file source.gene_enricher.hgnc_complete_set.json.dvc --yes \
  -d source/gene_enricher/version.txt \
  -o source/gene_enricher/hgnc_complete_set.json \
  "curl --verbose --progress-bar --ipv4 --connect-timeout 8 --max-time 120 --retry 128 --ftp-ssl --disable-epsv --ftp-pasv ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/json/hgnc_complete_set.json --output source/gene_enricher/hgnc_complete_set.json"



# commit DVC's record to github
git add .gitignore source.gene_enricher.hgnc_complete_set.json.dvc
git commit -m "add hgnc_complete_set to DVC"

# commit the data to DVC's remote
$ dvc push
Preparing to push data to s3://bmeg/dvc
[##############################] 100% Collecting information
[##############################] 100% source/gene_enricher/hgnc_complete_set.json

# determine md5 hash and view the remote repository

$ cat source.gene_enricher.hgnc_complete_set.json.dvc
cmd: curl --verbose --progress-bar --ipv4 --connect-timeout 8 --max-time 120 --retry
  128 --ftp-ssl --disable-epsv --ftp-pasv ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/json/hgnc_complete_set.json
  --output source/gene_enricher/hgnc_complete_set.json
deps:
- md5: 6808ca805661622ad65ae014a4b2a094
  path: source/gene_enricher/version.txt
md5: 1bada0b72fb7d4b523a0b9cbea3ef49f
outs:
- cache: true
  md5: cc3abedb40a6795c34b739a6df2b140e
  path: source/gene_enricher/hgnc_complete_set.json

$ cephos s3 ls s3://bmeg/dvc/cc/3abedb40a6795c34b739a6df2b140e
2018-12-18 21:50:34   31267685 3abedb40a6795c34b739a6df2b140e

# note that the file name ~ the md5 hash
```

Provenance
------

The dvc files were created using the following commands.
To recreate this section, run `python transform/dvc/dvc2cmd.py`

```

#
dvc run --file source.ccle.CCLE_DepMap_18q3_maf_20180718.txt.dvc --yes \
  -d source/ccle/version.txt \
  -o source/ccle/CCLE_DepMap_18q3_maf_20180718.txt \
  "wget https://data.broadinstitute.org/ccle/CCLE_DepMap_18q3_maf_20180718.txt -O source/ccle/CCLE_DepMap_18q3_maf_20180718.txt"
#
dvc run --file source.gene_enricher.hgnc_complete_set.json.dvc --yes \
  -d source/gene_enricher/version.txt \
  -o source/gene_enricher/hgnc_complete_set.json \
  "curl --verbose --progress-bar --ipv4 --connect-timeout 8 --max-time 120 --retry 128 --ftp-ssl --disable-epsv --ftp-pasv ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/json/hgnc_complete_set.json --output source/gene_enricher/hgnc_complete_set.json"
#
dvc run --file source.g2p.all.dvc --yes \
  -d source/g2p/version.txt \
  -o source/g2p/all.json \
  "curl https://s3-us-west-2.amazonaws.com/g2p-0.11/all.json | gzip > source/g2p/all.json.gz"
#
dvc import --file source.mc3.v0.2.8.PUBLIC.maf.gz.dvc --yes \
 remote://source/mc3/mc3.v0.2.8.PUBLIC.maf.gz\
 source/mc3/mc3.v0.2.8.PUBLIC.maf.gz
#
dvc import --file source.myvariants.biothings_current_old_hg19.json.gz.dvc --yes \
 remote://source/myvariant.info/biothings_current_old_hg19.json.gz\
 source/myvariant.info/biothings_current_old_hg19.json.gz
#
dvc import --file source.myvariant.info.harvested.json.gz.dvc --yes \
 remote://source/myvariant.info/harvested_myvariantinfo.json.gz\
 source/myvariant.info/harvested_myvariantinfo.json.gz
#
dvc import --file source.myvariant.info.metadata.fields.json.dvc --yes \
 remote://source/myvariant.info/metadata.fields.json\
 source/myvariant.info/metadata.fields.json
#
dvc import --file source.myvariant.info.metadata.json.dvc --yes \
 remote://source/myvariant.info/metadata.json\
 source/myvariant.info/metadata.json
#
dvc run --file source.ccle.DepMap-2018q3-celllines.csv.dvc --yes \
  -d source/ccle/version.txt \
  -o source/ccle/DepMap-2018q3-celllines.csv \
  "wget https://depmap.org/portal/download/api/download/external?file_name=processed_portal_downloads%2Fdepmap-public-2018q3-cell-line-metadata-0bd2.2%2FDepMap-2018q3-celllines.csv -O source/ccle/DepMap-2018q3-celllines.csv"
#
dvc run --file source.ccle.CCLE_DepMap_18q3_RNAseq_RPKM_20180718.gct.dvc --yes \
  -d source/ccle/version.txt \
  -o source/ccle/CCLE_DepMap_18q3_RNAseq_RPKM_20180718.gct \
  "wget https://data.broadinstitute.org/ccle/CCLE_DepMap_18q3_RNAseq_RPKM_20180718.gct -O source/ccle/CCLE_DepMap_18q3_RNAseq_RPKM_20180718.gct"
#
dvc run --file source.ccle.CCLE_NP24.2009_Drug_data_2015.02.24.csv.dvc --yes \
  -d source/ccle/version.txt \
  -o source/ccle/CCLE_NP24.2009_Drug_data_2015.02.24.csv \
  "wget https://data.broadinstitute.org/ccle_legacy_data/pharmacological_profiling/CCLE_NP24.2009_Drug_data_2015.02.24.csv -O source/ccle/CCLE_NP24.2009_Drug_data_2015.02.24.csv"
#
dvc run --file source.ccle.CCLE_tpm.tsv.gz.dvc --yes \
  -d source/ccle/version.txt \
  -o source/ccle/expression/CCLE_tpm.tsv.gz \
  "wget https://osf.io/brkh6/download -O source/ccle/expression/CCLE_tpm.tsv.gz"
#
dvc run --file source.gdsc.GDSC_AUC.csv.dvc --yes \
  -d source/gdsc/version.txt \
  -o source/gdsc/GDSC_AUC.csv \
  "wget https://depmap.org/portal/download/api/download/external?file_name=processed_portal_downloads%2Fsanger-gdsc-545e.1%2FGDSC_AUC.csv -O source/gdsc/GDSC_AUC.csv"
#
dvc run --file source.ensembl-protein.homo_sapiens.json.dvc --yes \
  -d source/ensembl-protein/version.txt \
  -o source/ensembl-protein/homo_sapiens.json \
  "transform/ensembl-protein/curl-homo_sapien.sh"
#
dvc run --file source.ensembl.Homo_sapiens.GRCh37.87.chr_patch_hapl_scaff.gff3.gz.dvc --yes \
  -d source/ensembl/version.txt \
  -o source/ensembl/Homo_sapiens.GRCh37.87.chr_patch_hapl_scaff.gff3.gz \
  "curl --verbose --progress-bar --ipv4 --connect-timeout 8 --max-time 120 --retry 128 --ftp-ssl --disable-epsv --ftp-pasv ftp://ftp.ensembl.org/pub/grch37/release-94/gff3/homo_sapiens/Homo_sapiens.GRCh37.87.chr_patch_hapl_scaff.gff3.gz --output source/ensembl/Homo_sapiens.GRCh37.87.chr_patch_hapl_scaff.gff3.gz"
#
dvc run --file source.goa_human.gaf.gz.dvc --yes \
  -d source/go/version.txt \
  -o source/go/goa_human.gaf.gz \
  "wget http://www.geneontology.org/gene-associations/goa_human.gaf.gz -O source/go/goa_human.gaf.gz"
#
dvc run --file source.go.HUMAN_9606_idmapping.dat.gz.dvc --yes \
  -d source/go/version.txt \
  -o source/go/HUMAN_9606_idmapping.dat.gz \
  "curl --verbose --progress-bar --ipv4 --connect-timeout 8 --max-time 120 --retry 128 --ftp-ssl --disable-epsv --ftp-pasv ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping.dat.gz --output source/go/HUMAN_9606_idmapping.dat.gz"
#
dvc run --file source.go.obo.dvc --yes \
  -d source/go/version.txt \
  -o source/go/go.obo \
  "wget http://purl.obolibrary.org/obo/go.obo -O source/go/go.obo"
#
dvc run --file source.gtex.GTEx_v7_Annotations_SampleAttributesDS.txt.dvc --yes \
  -d source/gtex/version.txt \
  -o source/gtex/GTEx_v7_Annotations_SampleAttributesDS.txt \
  "wget https://storage.googleapis.com/gtex_analysis_v7/annotations/GTEx_v7_Annotations_SampleAttributesDS.txt -O source/gtex/GTEx_v7_Annotations_SampleAttributesDS.txt"
#
dvc run --file source.gtex.GTEx_v7_Annotations_SubjectPhenotypesDS.txt.dvc --yes \
  -d source/gtex/version.txt \
  -o source/gtex/GTEx_v7_Annotations_SubjectPhenotypesDS.txt \
  "wget https://storage.googleapis.com/gtex_analysis_v7/annotations/GTEx_v7_Annotations_SubjectPhenotypesDS.txt -O source/gtex/GTEx_v7_Annotations_SubjectPhenotypesDS.txt"
#
dvc run --file source.gtex.GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct.dvc --yes \
  -d source/gtex/version.txt \
  -o source/gtex/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct \
  "wget https://storage.googleapis.com/gtex_analysis_v7/rna_seq_data/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct.gz -O source/gtex/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct"
#
dvc run --file source.pfam.clans.tsv.dvc --yes \
  -d source/pfam/version.txt \
  -o source/pfam/clans.tsv \
  "wget http://pfam.xfam.org/clans?output=text -O source/pfam/clans.tsv"
#
dvc run --file source.pfam.tar.gz.dvc --yes \
  -d transform/pfam/list.py \
  -o source/pfam/pfam.tar.gz \
  "python3 transform/pfam/download.py --archive ; tar xvfz source/pfam/pfam.tar.gz -C source/pfam"
#
dvc run --file source.pfam.id_list.txt.dvc --yes \
  -d source/pfam/version.txt \
  -o source/pfam/id_list.txt \
  "python3 transform/pfam/list.py -O source/pfam/id_list.txt"
#
dvc run --file source.tcga.TCGA_ID_MAP.csv.dvc --yes \
  -d source/tcga/version.txt \
  -o source/tcga/expression/transcript-level/TCGA_ID_MAP.csv \
  "wget https://osf.io/7qpsg/download -O source/tcga/expression/transcript-level/TCGA_ID_MAP.csv"
#
dvc run --file source.tcga.tcga-genomics.zip.dvc --yes \
  -d source/tcga/version.txt \
  -o source/tcga/expression/transcript-level/tcga-genomics.zip \
  "wget https://files.osf.io/v1/resources/gqrz9/providers/osfstorage/578518706c613b01f0a8325f/?zip= -O source/tcga/expression/transcript-level/tcga-genomics.zip"
#
dvc run --file source.tcga.TCGA_expression_tpm.tsv.gz.dvc --yes \
  -d source/tcga/expression/transcript-level/tcga-genomics.zip \
  -o source/tcga/expression/transcript-level/TCGA_ACC_tpm.tsv.gz \
  -o source/tcga/expression/transcript-level/TCGA_BLCA_tpm.tsv.gz \
  -o source/tcga/expression/transcript-level/TCGA_BRCA_tpm.tsv.gz \
  -o source/tcga/expression/transcript-level/TCGA_CESC_tpm.tsv.gz \
  -o source/tcga/expression/transcript-level/TCGA_CHOL_tpm.tsv.gz \
  -o source/tcga/expression/transcript-level/TCGA_COAD_tpm.tsv.gz \
  -o source/tcga/expression/transcript-level/TCGA_DLBC_tpm.tsv.gz \
  -o source/tcga/expression/transcript-level/TCGA_ESCA_tpm.tsv.gz \
  -o source/tcga/expression/transcript-level/TCGA_GBM_tpm.tsv.gz \
  -o source/tcga/expression/transcript-level/TCGA_HNSC_tpm.tsv.gz \
  -o source/tcga/expression/transcript-level/TCGA_KICH_tpm.tsv.gz \
  -o source/tcga/expression/transcript-level/TCGA_KIRC_tpm.tsv.gz \
  -o source/tcga/expression/transcript-level/TCGA_KIRP_tpm.tsv.gz \
  -o source/tcga/expression/transcript-level/TCGA_LAML_tpm.tsv.gz \
  -o source/tcga/expression/transcript-level/TCGA_LGG_tpm.tsv.gz \
  -o source/tcga/expression/transcript-level/TCGA_LIHC_tpm.tsv.gz \
  -o source/tcga/expression/transcript-level/TCGA_LUAD_tpm.tsv.gz \
  -o source/tcga/expression/transcript-level/TCGA_LUSC_tpm.tsv.gz \
  -o source/tcga/expression/transcript-level/TCGA_MESO_tpm.tsv.gz \
  -o source/tcga/expression/transcript-level/TCGA_OV_tpm.tsv.gz \
  -o source/tcga/expression/transcript-level/TCGA_PAAD_tpm.tsv.gz \
  -o source/tcga/expression/transcript-level/TCGA_PCPG_tpm.tsv.gz \
  -o source/tcga/expression/transcript-level/TCGA_PRAD_tpm.tsv.gz \
  -o source/tcga/expression/transcript-level/TCGA_READ_tpm.tsv.gz \
  -o source/tcga/expression/transcript-level/TCGA_SARC_tpm.tsv.gz \
  -o source/tcga/expression/transcript-level/TCGA_SKCM_tpm.tsv.gz \
  -o source/tcga/expression/transcript-level/TCGA_STAD_tpm.tsv.gz \
  -o source/tcga/expression/transcript-level/TCGA_TGCT_tpm.tsv.gz \
  -o source/tcga/expression/transcript-level/TCGA_THCA_tpm.tsv.gz \
  -o source/tcga/expression/transcript-level/TCGA_THYM_tpm.tsv.gz \
  -o source/tcga/expression/transcript-level/TCGA_UCEC_tpm.tsv.gz \
  -o source/tcga/expression/transcript-level/TCGA_UCS_tpm.tsv.gz \
  -o source/tcga/expression/transcript-level/TCGA_UVM_tpm.tsv.gz \
  "unzip -j source/tcga/expression/transcript-level/tcga-genomics.zip *_tpm.tsv.gz -d source/tcga/expression/transcript-level"
#
dvc run --file outputs.ccle.maf.dvc --yes \
  -d source/ccle/CCLE_DepMap_18q3_maf_20180718.txt \
  -d source/gene_enricher/hgnc_complete_set.json \
  -o outputs/ccle/Allele.Vertex.json.gz \
  -o outputs/ccle/AlleleCall.Edge.json.gz \
  -o outputs/ccle/AlleleIn.Edge.json.gz \
  -o outputs/ccle/Callset.Vertex.json.gz \
  -o outputs/ccle/CallsetFor.Edge.json.gz \
  "python3 transform/ccle/ccle_maf_transform.py"
#
dvc run --file outputs.g2p.dvc --yes \
  -d source/g2p/all.json \
  -d source/gene_enricher/hgnc_complete_set.json \
  -o outputs/g2p/Allele.Vertex.json.gz \
  -o outputs/g2p/AlleleIn.Edge.json.gz \
  -o outputs/g2p/Compound.Vertex.json.gz \
  -o outputs/g2p/Deadletter.Vertex.json.gz \
  -o outputs/g2p/G2PAssociation.Vertex.json.gz \
  -o outputs/g2p/HasAlleleFeature.Edge.json.gz \
  -o outputs/g2p/HasEnvironment.Edge.json.gz \
  -o outputs/g2p/HasGeneFeature.Edge.json.gz \
  -o outputs/g2p/HasMinimalAlleleFeature.Edge.json.gz \
  -o outputs/g2p/HasPhenotype.Edge.json.gz \
  -o outputs/g2p/HasSupportingReference.Edge.json.gz \
  -o outputs/g2p/MinimalAllele.Vertex.json.gz \
  -o outputs/g2p/MinimalAlleleIn.Edge.json.gz \
  -o outputs/g2p/Phenotype.Vertex.json.gz \
  "python3 -m transform.g2p.transform"
#
dvc run --file outputs.gdc.cases.dvc --yes \
  -d source/gdc/version.txt \
  -o outputs/gdc/Aliquot.Vertex.json.gz \
  -o outputs/gdc/AliquotFor.Edge.json.gz \
  -o outputs/gdc/Biosample.Vertex.json.gz \
  -o outputs/gdc/BiosampleFor.Edge.json.gz \
  -o outputs/gdc/Compound.Vertex.json.gz \
  -o outputs/gdc/InProject.Edge.json.gz \
  -o outputs/gdc/Individual.Vertex.json.gz \
  -o outputs/gdc/TreatedWith.Edge.json.gz \
  "python3 -m transform.gdc.cases"
#
dvc run --file outputs.mc3.dvc --yes \
  -d outputs/gdc/Aliquot.Vertex.json.gz \
  -d source/mc3/mc3.v0.2.8.PUBLIC.maf.gz \
  -o outputs/mc3/Allele.Vertex.json.gz \
  -o outputs/mc3/AlleleCall.Edge.json.gz \
  -o outputs/mc3/AlleleIn.Edge.json.gz \
  -o outputs/mc3/Callset.Vertex.json.gz \
  -o outputs/mc3/CallsetFor.Edge.json.gz \
  -o outputs/mc3/Deadletter.Vertex.json.gz \
  "python3 transform/mc3/mc3_maf_transform.py"
#
dvc run --file outputs.allele.dvc --yes \
  -d outputs/ccle/Allele.Vertex.json.gz \
  -d outputs/g2p/Allele.Vertex.json.gz \
  -d outputs/g2p/MinimalAllele.Vertex.json.gz \
  -d outputs/mc3/Allele.Vertex.json.gz \
  -d outputs/myvariant.info/myvariant.info.Allele.Vertex.json.gz \
  -d source/myvariant.info/biothings_current_old_hg19.json.gz \
  -d source/myvariant.info/harvested_myvariantinfo.json.gz \
  -d source/myvariant.info/metadata.fields.json \
  -d source/myvariant.info/metadata.json \
  -o outputs/allele/Allele.Vertex.json.gz \
  "python3 transform/allele/transform.py"
#
dvc run --file outputs.ccle.dvc --yes \
  -d source/ccle/DepMap-2018q3-celllines.csv \
  -o outputs/ccle/Aliquot.Vertex.json.gz \
  -o outputs/ccle/AliquotFor.Edge.json.gz \
  -o outputs/ccle/Biosample.Vertex.json.gz \
  -o outputs/ccle/BiosampleFor.Edge.json.gz \
  -o outputs/ccle/InProject.Edge.json.gz \
  -o outputs/ccle/Individual.Vertex.json.gz \
  -o outputs/ccle/Phenotype.Vertex.json.gz \
  -o outputs/ccle/PhenotypeOf.Edge.json.gz \
  -o outputs/ccle/Project.Vertex.json.gz \
  "python3 transform/ccle/samples.py"
#
dvc run --file outputs.ccle.expression.dvc --yes \
  -d source/ccle/CCLE_DepMap_18q3_RNAseq_RPKM_20180718.gct \
  -o outputs/ccle/Expression.Vertex.json.gz \
  -o outputs/ccle/ExpressionOf.Edge.json.gz \
  "python3 transform/ccle/expression.py"
#
dvc run --file outputs.ccle.drug_response.dvc --yes \
  -d source/ccle/CCLE_NP24.2009_Drug_data_2015.02.24.csv \
  -d outputs/ccle/Biosample.Vertex.json.gz \
  -o outputs/ccle/drug_response.Aliquot.Vertex.json.gz \
  -o outputs/ccle/drug_response.AliquotFor.Edge.json.gz \
  -o outputs/ccle/drug_response.Biosample.Vertex.json.gz \
  -o outputs/ccle/drug_response.BiosampleFor.Edge.json.gz \
  -o outputs/ccle/drug_response.Compound.Vertex.json.gz \
  -o outputs/ccle/drug_response.InProject.Edge.json.gz \
  -o outputs/ccle/drug_response.Individual.Vertex.json.gz \
  -o outputs/ccle/drug_response.ParamacalogicalProfile.Vertex.json.gz \
  -o outputs/ccle/drug_response.ParamacalogicalProfileIn.Edge.json.gz \
  -o outputs/ccle/drug_response.Project.Vertex.json.gz \
  -o outputs/ccle/drug_response.ResponseTo.Edge.json.gz \
  "python3 transform/ccle/drug_response.py"
#
dvc run --file outputs.ccle.expression_tatlow.dvc --yes \
  -d source/ccle/expression/CCLE_tpm.tsv.gz \
  -d outputs/ccle/Biosample.Vertex.json.gz \
  -d source/tcga/expression/transcript-level/TCGA_ID_MAP.csv \
  -o outputs/ccle/tatlow.Aliquot.Vertex.json.gz \
  -o outputs/ccle/tatlow.AliquotFor.Edge.json.gz \
  -o outputs/ccle/tatlow.Biosample.Vertex.json.gz \
  -o outputs/ccle/tatlow.BiosampleFor.Edge.json.gz \
  -o outputs/ccle/tatlow.Expression.Vertex.json.gz \
  -o outputs/ccle/tatlow.ExpressionOf.Edge.json.gz \
  -o outputs/ccle/tatlow.InProject.Edge.json.gz \
  -o outputs/ccle/tatlow.Individual.Vertex.json.gz \
  -o outputs/ccle/tatlow.Project.Vertex.json.gz \
  "python3 transform/ccle/expression_tatlow.py"
#
dvc run --file outputs.gdsc.dvc --yes \
  -d source/gdsc/GDSC_AUC.csv \
  -o outputs/gdsc/gdsc.Compound.Vertex.json.gz \
  -o outputs/gdsc/gdsc.DrugResponse.Vertex.json.gz \
  -o outputs/gdsc/gdsc.DrugResponseIn.Edge.json.gz \
  -o outputs/gdsc/gdsc.ResponseTo.Edge.json.gz \
  "python3 transform/gdsc/response.py"
#
dvc run --file outputs.compound.normalized.dvc --yes \
  -d source/compound/sqlite.db \
  -d outputs/ccle/Compound.Vertex.json.gz \
  -d outputs/ccle/drug_response.Compound.Vertex.json.gz \
  -d outputs/g2p/Compound.Vertex.json.gz \
  -d outputs/gdc/Compound.Vertex.json.gz \
  -d outputs/gdsc/gdsc.Compound.Vertex.json.gz \
  -d outputs/g2p/HasEnvironment.Edge.json.gz \
  -d outputs/ccle/ResponseTo.Edge.json.gz \
  -d outputs/ccle/drug_response.ResponseTo.Edge.json.gz \
  -d outputs/gdsc/gdsc.ResponseTo.Edge.json.gz \
  -d outputs/gdc/TreatedWith.Edge.json.gz \
  -o outputs/compound/normalized.Compound.Vertex.json.gz \
  -o outputs/compound/normalized.HasEnvironment.Edge.json.gz \
  -o outputs/compound/normalized.ResponseTo.Edge.json.gz \
  -o outputs/compound/normalized.TreatedWith.Edge.json.gz \
  "python3 transform/compound/transform.py"
#
dvc run --file outputs.ensembl-protein.dvc --yes \
  -d source/ensembl-protein/homo_sapiens.json \
  -o outputs/ensembl-protein/PFAMAlignment.Edge.json.gz \
  -o outputs/ensembl-protein/Protein.Vertex.json.gz \
  -o outputs/ensembl-protein/ProteinFor.Edge.json.gz \
  -o outputs/ensembl-protein/ProteinStructure.Vertex.json.gz \
  -o outputs/ensembl-protein/StructureFor.Edge.json.gz \
  "python3 transform/ensembl-protein/ensembl-protein-transform.py"
#
dvc run --file outputs.ensembl.dvc --yes \
  -d source/ensembl/Homo_sapiens.GRCh37.87.chr_patch_hapl_scaff.gff3.gz \
  -o outputs/ensembl/Exon.Vertex.json.gz \
  -o outputs/ensembl/ExonFor.Edge.json.gz \
  -o outputs/ensembl/Gene.Vertex.json.gz \
  -o outputs/ensembl/Transcript.Vertex.json.gz \
  -o outputs/ensembl/TranscriptFor.Edge.json.gz \
  "python3 transform/ensembl/transform.py"
#
dvc run --file outputs.ensembl.missing_transcripts.dvc --yes \
  -d source/ensembl/missing_transcript_ids.txt \
  -o outputs/ensembl/missing.Transcript.Vertex.json.gz \
  "python3 transform/ensembl/missing_transcripts.py"
#
dvc run --file outputs.gdc.projects.dvc --yes \
  -d source/gdc/version.txt \
  -o outputs/gdc/Project.Vertex.json.gz \
  "python3 transform/gdc/projects.py"
#
dvc run --file outputs.go.gaf2schema.dvc --yes \
  -d source/go/goa_human.gaf.gz \
  -d source/go/HUMAN_9606_idmapping.dat.gz \
  -o outputs/go/GeneOntologyAnnotation.Edge.json.gz \
  "python3 transform/go/go_gaf2schema.py source/go/goa_human.gaf.gz source/go/HUMAN_9606_idmapping.dat.gz go"
#
dvc run --file outputs.go.dvc --yes \
  -d source/go/go.obo \
  -o outputs/go/GeneOntologyIsA.Edge.json.gz \
  -o outputs/go/GeneOntologyTerm.Vertex.json.gz \
  "python3 transform/go/go_obo2schema.py source/go/go.obo go"
#
dvc run --file outputs.gtex.dvc --yes \
  -d source/gtex/GTEx_v7_Annotations_SampleAttributesDS.txt \
  -d source/gtex/GTEx_v7_Annotations_SubjectPhenotypesDS.txt \
  -o outputs/gtex/gtex.Aliquot.Vertex.json.gz \
  -o outputs/gtex/gtex.AliquotFor.Edge.json.gz \
  -o outputs/gtex/gtex.Biosample.Vertex.json.gz \
  -o outputs/gtex/gtex.BiosampleFor.Edge.json.gz \
  -o outputs/gtex/gtex.InProject.Edge.json.gz \
  -o outputs/gtex/gtex.Individual.Vertex.json.gz \
  -o outputs/gtex/gtex.Project.Vertex.json.gz \
  "python3 transform/gtex/samples.py"
#
dvc run --file outputs.gtex.expression.dvc --yes \
  -d source/gtex/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct \
  -o outputs/gtex/gtex.Expression.Vertex.json.gz \
  -o outputs/gtex/gtex.ExpressionOf.Edge.json.gz \
  "python3 transform/gtex/expression.py"
#
dvc run --file outputs.pfam.dvc --yes \
  -d source/pfam/clans.tsv \
  -d source/pfam/pfam.tar.gz \
  -d source/pfam/id_list.txt \
  -o outputs/pfam/GeneOntologyAnnotation.Edge.json.gz \
  -o outputs/pfam/PFAMClan.Vertex.json.gz \
  -o outputs/pfam/PFAMClanMember.Edge.json.gz \
  -o outputs/pfam/PFAMFamily.Vertex.json.gz \
  "python3 transform/pfam/transform.py"
#
dvc run --file outputs.phenotype.dvc --yes \
  -d outputs/ccle/Phenotype.Vertex.json.gz \
  -d outputs/g2p/Phenotype.Vertex.json.gz \
  -o outputs/phenotype/normalized.HasPhenotype.Edge.json.gz \
  -o outputs/phenotype/normalized.Phenotype.Vertex.json.gz \
  -o outputs/phenotype/normalized.PhenotypeOf.Edge.json.gz \
  "python3 transform/phenotype/transform.py"
#
dvc run --file outputs.publication.dvc --yes \
  -d outputs/g2p/HasSupportingReference.Edge.json.gz \
  -o outputs/publication/stub.Publication.Vertex.json.gz \
  "python3 transform/publication/transform.py"
#
dvc run --file outputs.tcga.expression.dvc --yes \
  -d source/tcga/expression/transcript-level/TCGA_ID_MAP.csv \
  -d source/tcga/expression/transcript-level/TCGA_ACC_tpm.tsv.gz \
  -d source/tcga/expression/transcript-level/TCGA_BLCA_tpm.tsv.gz \
  -d source/tcga/expression/transcript-level/TCGA_BRCA_tpm.tsv.gz \
  -d source/tcga/expression/transcript-level/TCGA_CESC_tpm.tsv.gz \
  -d source/tcga/expression/transcript-level/TCGA_CHOL_tpm.tsv.gz \
  -d source/tcga/expression/transcript-level/TCGA_COAD_tpm.tsv.gz \
  -d source/tcga/expression/transcript-level/TCGA_DLBC_tpm.tsv.gz \
  -d source/tcga/expression/transcript-level/TCGA_ESCA_tpm.tsv.gz \
  -d source/tcga/expression/transcript-level/TCGA_GBM_tpm.tsv.gz \
  -d source/tcga/expression/transcript-level/TCGA_HNSC_tpm.tsv.gz \
  -d source/tcga/expression/transcript-level/TCGA_KICH_tpm.tsv.gz \
  -d source/tcga/expression/transcript-level/TCGA_KIRC_tpm.tsv.gz \
  -d source/tcga/expression/transcript-level/TCGA_KIRP_tpm.tsv.gz \
  -d source/tcga/expression/transcript-level/TCGA_LAML_tpm.tsv.gz \
  -d source/tcga/expression/transcript-level/TCGA_LGG_tpm.tsv.gz \
  -d source/tcga/expression/transcript-level/TCGA_LIHC_tpm.tsv.gz \
  -d source/tcga/expression/transcript-level/TCGA_LUAD_tpm.tsv.gz \
  -d source/tcga/expression/transcript-level/TCGA_LUSC_tpm.tsv.gz \
  -d source/tcga/expression/transcript-level/TCGA_MESO_tpm.tsv.gz \
  -d source/tcga/expression/transcript-level/TCGA_OV_tpm.tsv.gz \
  -d source/tcga/expression/transcript-level/TCGA_PAAD_tpm.tsv.gz \
  -d source/tcga/expression/transcript-level/TCGA_PCPG_tpm.tsv.gz \
  -d source/tcga/expression/transcript-level/TCGA_PRAD_tpm.tsv.gz \
  -d source/tcga/expression/transcript-level/TCGA_READ_tpm.tsv.gz \
  -d source/tcga/expression/transcript-level/TCGA_SARC_tpm.tsv.gz \
  -d source/tcga/expression/transcript-level/TCGA_SKCM_tpm.tsv.gz \
  -d source/tcga/expression/transcript-level/TCGA_STAD_tpm.tsv.gz \
  -d source/tcga/expression/transcript-level/TCGA_TGCT_tpm.tsv.gz \
  -d source/tcga/expression/transcript-level/TCGA_THCA_tpm.tsv.gz \
  -d source/tcga/expression/transcript-level/TCGA_THYM_tpm.tsv.gz \
  -d source/tcga/expression/transcript-level/TCGA_UCEC_tpm.tsv.gz \
  -d source/tcga/expression/transcript-level/TCGA_UCS_tpm.tsv.gz \
  -d source/tcga/expression/transcript-level/TCGA_UVM_tpm.tsv.gz \
  -o outputs/tcga/ACC.Expression.Vertex.json.gz \
  -o outputs/tcga/ACC.ExpressionOf.Edge.json.gz \
  -o outputs/tcga/BLCA.Expression.Vertex.json.gz \
  -o outputs/tcga/BLCA.ExpressionOf.Edge.json.gz \
  -o outputs/tcga/BRCA.Expression.Vertex.json.gz \
  -o outputs/tcga/BRCA.ExpressionOf.Edge.json.gz \
  -o outputs/tcga/CESC.Expression.Vertex.json.gz \
  -o outputs/tcga/CESC.ExpressionOf.Edge.json.gz \
  -o outputs/tcga/CHOL.Expression.Vertex.json.gz \
  -o outputs/tcga/CHOL.ExpressionOf.Edge.json.gz \
  -o outputs/tcga/COAD.Expression.Vertex.json.gz \
  -o outputs/tcga/COAD.ExpressionOf.Edge.json.gz \
  -o outputs/tcga/DLBC.Expression.Vertex.json.gz \
  -o outputs/tcga/DLBC.ExpressionOf.Edge.json.gz \
  -o outputs/tcga/ESCA.Expression.Vertex.json.gz \
  -o outputs/tcga/ESCA.ExpressionOf.Edge.json.gz \
  -o outputs/tcga/GBM.Expression.Vertex.json.gz \
  -o outputs/tcga/GBM.ExpressionOf.Edge.json.gz \
  -o outputs/tcga/HNSC.Expression.Vertex.json.gz \
  -o outputs/tcga/HNSC.ExpressionOf.Edge.json.gz \
  -o outputs/tcga/KICH.Expression.Vertex.json.gz \
  -o outputs/tcga/KICH.ExpressionOf.Edge.json.gz \
  -o outputs/tcga/KIRC.Expression.Vertex.json.gz \
  -o outputs/tcga/KIRC.ExpressionOf.Edge.json.gz \
  -o outputs/tcga/KIRP.Expression.Vertex.json.gz \
  -o outputs/tcga/KIRP.ExpressionOf.Edge.json.gz \
  -o outputs/tcga/LAML.Expression.Vertex.json.gz \
  -o outputs/tcga/LAML.ExpressionOf.Edge.json.gz \
  -o outputs/tcga/LGG.Expression.Vertex.json.gz \
  -o outputs/tcga/LGG.ExpressionOf.Edge.json.gz \
  -o outputs/tcga/LIHC.Expression.Vertex.json.gz \
  -o outputs/tcga/LIHC.ExpressionOf.Edge.json.gz \
  -o outputs/tcga/LUAD.Expression.Vertex.json.gz \
  -o outputs/tcga/LUAD.ExpressionOf.Edge.json.gz \
  -o outputs/tcga/LUSC.Expression.Vertex.json.gz \
  -o outputs/tcga/LUSC.ExpressionOf.Edge.json.gz \
  -o outputs/tcga/MESO.Expression.Vertex.json.gz \
  -o outputs/tcga/MESO.ExpressionOf.Edge.json.gz \
  -o outputs/tcga/OV.Expression.Vertex.json.gz \
  -o outputs/tcga/OV.ExpressionOf.Edge.json.gz \
  -o outputs/tcga/PAAD.Expression.Vertex.json.gz \
  -o outputs/tcga/PAAD.ExpressionOf.Edge.json.gz \
  -o outputs/tcga/PCPG.Expression.Vertex.json.gz \
  -o outputs/tcga/PCPG.ExpressionOf.Edge.json.gz \
  -o outputs/tcga/PRAD.Expression.Vertex.json.gz \
  -o outputs/tcga/PRAD.ExpressionOf.Edge.json.gz \
  -o outputs/tcga/READ.Expression.Vertex.json.gz \
  -o outputs/tcga/READ.ExpressionOf.Edge.json.gz \
  -o outputs/tcga/SARC.Expression.Vertex.json.gz \
  -o outputs/tcga/SARC.ExpressionOf.Edge.json.gz \
  -o outputs/tcga/SKCM.Expression.Vertex.json.gz \
  -o outputs/tcga/SKCM.ExpressionOf.Edge.json.gz \
  -o outputs/tcga/STAD.Expression.Vertex.json.gz \
  -o outputs/tcga/STAD.ExpressionOf.Edge.json.gz \
  -o outputs/tcga/TGCT.Expression.Vertex.json.gz \
  -o outputs/tcga/TGCT.ExpressionOf.Edge.json.gz \
  -o outputs/tcga/THCA.Expression.Vertex.json.gz \
  -o outputs/tcga/THCA.ExpressionOf.Edge.json.gz \
  -o outputs/tcga/THYM.Expression.Vertex.json.gz \
  -o outputs/tcga/THYM.ExpressionOf.Edge.json.gz \
  -o outputs/tcga/UCEC.Expression.Vertex.json.gz \
  -o outputs/tcga/UCEC.ExpressionOf.Edge.json.gz \
  -o outputs/tcga/UCS.Expression.Vertex.json.gz \
  -o outputs/tcga/UCS.ExpressionOf.Edge.json.gz \
  -o outputs/tcga/UVM.Expression.Vertex.json.gz \
  -o outputs/tcga/UVM.ExpressionOf.Edge.json.gz \
  "python3 transform/tcga/expression.py"
#
dvc run --file outputs.bmeg_manifest.dvc --yes \
  -d outputs/allele/Allele.Vertex.json.gz \
  -d outputs/ccle/Aliquot.Vertex.json.gz \
  -d outputs/ccle/AliquotFor.Edge.json.gz \
  -d outputs/ccle/AlleleCall.Edge.json.gz \
  -d outputs/ccle/AlleleIn.Edge.json.gz \
  -d outputs/ccle/Biosample.Vertex.json.gz \
  -d outputs/ccle/BiosampleFor.Edge.json.gz \
  -d outputs/ccle/Callset.Vertex.json.gz \
  -d outputs/ccle/CallsetFor.Edge.json.gz \
  -d outputs/ccle/Expression.Vertex.json.gz \
  -d outputs/ccle/ExpressionOf.Edge.json.gz \
  -d outputs/ccle/InProject.Edge.json.gz \
  -d outputs/ccle/Individual.Vertex.json.gz \
  -d outputs/ccle/Project.Vertex.json.gz \
  -d outputs/ccle/drug_response.Aliquot.Vertex.json.gz \
  -d outputs/ccle/drug_response.AliquotFor.Edge.json.gz \
  -d outputs/ccle/drug_response.Biosample.Vertex.json.gz \
  -d outputs/ccle/drug_response.BiosampleFor.Edge.json.gz \
  -d outputs/ccle/drug_response.InProject.Edge.json.gz \
  -d outputs/ccle/drug_response.Individual.Vertex.json.gz \
  -d outputs/ccle/drug_response.ParamacalogicalProfile.Vertex.json.gz \
  -d outputs/ccle/drug_response.ParamacalogicalProfileIn.Edge.json.gz \
  -d outputs/ccle/drug_response.Project.Vertex.json.gz \
  -d outputs/ccle/tatlow.Aliquot.Vertex.json.gz \
  -d outputs/ccle/tatlow.AliquotFor.Edge.json.gz \
  -d outputs/ccle/tatlow.Biosample.Vertex.json.gz \
  -d outputs/ccle/tatlow.BiosampleFor.Edge.json.gz \
  -d outputs/ccle/tatlow.Expression.Vertex.json.gz \
  -d outputs/ccle/tatlow.ExpressionOf.Edge.json.gz \
  -d outputs/ccle/tatlow.InProject.Edge.json.gz \
  -d outputs/ccle/tatlow.Individual.Vertex.json.gz \
  -d outputs/ccle/tatlow.Project.Vertex.json.gz \
  -d outputs/compound/normalized.Compound.Vertex.json.gz \
  -d outputs/compound/normalized.HasEnvironment.Edge.json.gz \
  -d outputs/compound/normalized.ResponseTo.Edge.json.gz \
  -d outputs/compound/normalized.TreatedWith.Edge.json.gz \
  -d outputs/ensembl-protein/PFAMAlignment.Edge.json.gz \
  -d outputs/ensembl-protein/Protein.Vertex.json.gz \
  -d outputs/ensembl-protein/ProteinFor.Edge.json.gz \
  -d outputs/ensembl-protein/ProteinStructure.Vertex.json.gz \
  -d outputs/ensembl-protein/StructureFor.Edge.json.gz \
  -d outputs/ensembl/Exon.Vertex.json.gz \
  -d outputs/ensembl/ExonFor.Edge.json.gz \
  -d outputs/ensembl/Gene.Vertex.json.gz \
  -d outputs/ensembl/Transcript.Vertex.json.gz \
  -d outputs/ensembl/TranscriptFor.Edge.json.gz \
  -d outputs/ensembl/missing.Transcript.Vertex.json.gz \
  -d outputs/g2p/AlleleIn.Edge.json.gz \
  -d outputs/g2p/G2PAssociation.Vertex.json.gz \
  -d outputs/g2p/HasAlleleFeature.Edge.json.gz \
  -d outputs/g2p/HasGeneFeature.Edge.json.gz \
  -d outputs/g2p/HasMinimalAlleleFeature.Edge.json.gz \
  -d outputs/g2p/HasSupportingReference.Edge.json.gz \
  -d outputs/g2p/MinimalAllele.Vertex.json.gz \
  -d outputs/g2p/MinimalAlleleIn.Edge.json.gz \
  -d outputs/gdc/Aliquot.Vertex.json.gz \
  -d outputs/gdc/AliquotFor.Edge.json.gz \
  -d outputs/gdc/Biosample.Vertex.json.gz \
  -d outputs/gdc/BiosampleFor.Edge.json.gz \
  -d outputs/gdc/InProject.Edge.json.gz \
  -d outputs/gdc/Individual.Vertex.json.gz \
  -d outputs/gdc/Project.Vertex.json.gz \
  -d outputs/gdsc/gdsc.DrugResponse.Vertex.json.gz \
  -d outputs/gdsc/gdsc.DrugResponseIn.Edge.json.gz \
  -d outputs/go/GeneOntologyAnnotation.Edge.json.gz \
  -d outputs/go/GeneOntologyIsA.Edge.json.gz \
  -d outputs/go/GeneOntologyTerm.Vertex.json.gz \
  -d outputs/gtex/gtex.Aliquot.Vertex.json.gz \
  -d outputs/gtex/gtex.AliquotFor.Edge.json.gz \
  -d outputs/gtex/gtex.Biosample.Vertex.json.gz \
  -d outputs/gtex/gtex.BiosampleFor.Edge.json.gz \
  -d outputs/gtex/gtex.Expression.Vertex.json.gz \
  -d outputs/gtex/gtex.ExpressionOf.Edge.json.gz \
  -d outputs/gtex/gtex.InProject.Edge.json.gz \
  -d outputs/gtex/gtex.Individual.Vertex.json.gz \
  -d outputs/gtex/gtex.Project.Vertex.json.gz \
  -d outputs/mc3/AlleleCall.Edge.json.gz \
  -d outputs/mc3/AlleleIn.Edge.json.gz \
  -d outputs/mc3/Callset.Vertex.json.gz \
  -d outputs/mc3/CallsetFor.Edge.json.gz \
  -d outputs/pfam/GeneOntologyAnnotation.Edge.json.gz \
  -d outputs/pfam/PFAMClan.Vertex.json.gz \
  -d outputs/pfam/PFAMClanMember.Edge.json.gz \
  -d outputs/pfam/PFAMFamily.Vertex.json.gz \
  -d outputs/phenotype/normalized.HasPhenotype.Edge.json.gz \
  -d outputs/phenotype/normalized.Phenotype.Vertex.json.gz \
  -d outputs/phenotype/normalized.PhenotypeOf.Edge.json.gz \
  -d outputs/publication/stub.Publication.Vertex.json.gz \
  -d outputs/tcga/ACC.Expression.Vertex.json.gz \
  -d outputs/tcga/ACC.ExpressionOf.Edge.json.gz \
  -d outputs/tcga/BLCA.Expression.Vertex.json.gz \
  -d outputs/tcga/BLCA.ExpressionOf.Edge.json.gz \
  -d outputs/tcga/BRCA.Expression.Vertex.json.gz \
  -d outputs/tcga/BRCA.ExpressionOf.Edge.json.gz \
  -d outputs/tcga/CESC.Expression.Vertex.json.gz \
  -d outputs/tcga/CESC.ExpressionOf.Edge.json.gz \
  -d outputs/tcga/CHOL.Expression.Vertex.json.gz \
  -d outputs/tcga/CHOL.ExpressionOf.Edge.json.gz \
  -d outputs/tcga/COAD.Expression.Vertex.json.gz \
  -d outputs/tcga/COAD.ExpressionOf.Edge.json.gz \
  -d outputs/tcga/DLBC.Expression.Vertex.json.gz \
  -d outputs/tcga/DLBC.ExpressionOf.Edge.json.gz \
  -d outputs/tcga/ESCA.Expression.Vertex.json.gz \
  -d outputs/tcga/ESCA.ExpressionOf.Edge.json.gz \
  -d outputs/tcga/GBM.Expression.Vertex.json.gz \
  -d outputs/tcga/GBM.ExpressionOf.Edge.json.gz \
  -d outputs/tcga/HNSC.Expression.Vertex.json.gz \
  -d outputs/tcga/HNSC.ExpressionOf.Edge.json.gz \
  -d outputs/tcga/KICH.Expression.Vertex.json.gz \
  -d outputs/tcga/KICH.ExpressionOf.Edge.json.gz \
  -d outputs/tcga/KIRC.Expression.Vertex.json.gz \
  -d outputs/tcga/KIRC.ExpressionOf.Edge.json.gz \
  -d outputs/tcga/KIRP.Expression.Vertex.json.gz \
  -d outputs/tcga/KIRP.ExpressionOf.Edge.json.gz \
  -d outputs/tcga/LAML.Expression.Vertex.json.gz \
  -d outputs/tcga/LAML.ExpressionOf.Edge.json.gz \
  -d outputs/tcga/LGG.Expression.Vertex.json.gz \
  -d outputs/tcga/LGG.ExpressionOf.Edge.json.gz \
  -d outputs/tcga/LIHC.Expression.Vertex.json.gz \
  -d outputs/tcga/LIHC.ExpressionOf.Edge.json.gz \
  -d outputs/tcga/LUAD.Expression.Vertex.json.gz \
  -d outputs/tcga/LUAD.ExpressionOf.Edge.json.gz \
  -d outputs/tcga/LUSC.Expression.Vertex.json.gz \
  -d outputs/tcga/LUSC.ExpressionOf.Edge.json.gz \
  -d outputs/tcga/MESO.Expression.Vertex.json.gz \
  -d outputs/tcga/MESO.ExpressionOf.Edge.json.gz \
  -d outputs/tcga/OV.Expression.Vertex.json.gz \
  -d outputs/tcga/OV.ExpressionOf.Edge.json.gz \
  -d outputs/tcga/PAAD.Expression.Vertex.json.gz \
  -d outputs/tcga/PAAD.ExpressionOf.Edge.json.gz \
  -d outputs/tcga/PCPG.Expression.Vertex.json.gz \
  -d outputs/tcga/PCPG.ExpressionOf.Edge.json.gz \
  -d outputs/tcga/PRAD.Expression.Vertex.json.gz \
  -d outputs/tcga/PRAD.ExpressionOf.Edge.json.gz \
  -d outputs/tcga/READ.Expression.Vertex.json.gz \
  -d outputs/tcga/READ.ExpressionOf.Edge.json.gz \
  -d outputs/tcga/SARC.Expression.Vertex.json.gz \
  -d outputs/tcga/SARC.ExpressionOf.Edge.json.gz \
  -d outputs/tcga/SKCM.Expression.Vertex.json.gz \
  -d outputs/tcga/SKCM.ExpressionOf.Edge.json.gz \
  -d outputs/tcga/STAD.Expression.Vertex.json.gz \
  -d outputs/tcga/STAD.ExpressionOf.Edge.json.gz \
  -d outputs/tcga/TGCT.Expression.Vertex.json.gz \
  -d outputs/tcga/TGCT.ExpressionOf.Edge.json.gz \
  -d outputs/tcga/THCA.Expression.Vertex.json.gz \
  -d outputs/tcga/THCA.ExpressionOf.Edge.json.gz \
  -d outputs/tcga/THYM.Expression.Vertex.json.gz \
  -d outputs/tcga/THYM.ExpressionOf.Edge.json.gz \
  -d outputs/tcga/UCEC.Expression.Vertex.json.gz \
  -d outputs/tcga/UCEC.ExpressionOf.Edge.json.gz \
  -d outputs/tcga/UCS.Expression.Vertex.json.gz \
  -d outputs/tcga/UCS.ExpressionOf.Edge.json.gz \
  -d outputs/tcga/UVM.Expression.Vertex.json.gz \
  -d outputs/tcga/UVM.ExpressionOf.Edge.json.gz \
  -d outputs/tcga/TCGA-ACC.CopyNumberAlteration.Vertex.json.gz \
  -d outputs/tcga/TCGA-ACC.CopyNumberAlterationOf.Edge.json.gz \
  -d outputs/tcga/TCGA-BLCA.CopyNumberAlteration.Vertex.json.gz \
  -d outputs/tcga/TCGA-BLCA.CopyNumberAlterationOf.Edge.json.gz \
  -d outputs/tcga/TCGA-BRCA.CopyNumberAlteration.Vertex.json.gz \
  -d outputs/tcga/TCGA-BRCA.CopyNumberAlterationOf.Edge.json.gz \
  -d outputs/tcga/TCGA-CESC.CopyNumberAlteration.Vertex.json.gz \
  -d outputs/tcga/TCGA-CESC.CopyNumberAlterationOf.Edge.json.gz \
  -d outputs/tcga/TCGA-CHOL.CopyNumberAlteration.Vertex.json.gz \
  -d outputs/tcga/TCGA-CHOL.CopyNumberAlterationOf.Edge.json.gz \
  -d outputs/tcga/TCGA-COAD.CopyNumberAlteration.Vertex.json.gz \
  -d outputs/tcga/TCGA-COAD.CopyNumberAlterationOf.Edge.json.gz \
  -d outputs/tcga/TCGA-DLBC.CopyNumberAlteration.Vertex.json.gz \
  -d outputs/tcga/TCGA-DLBC.CopyNumberAlterationOf.Edge.json.gz \
  -d outputs/tcga/TCGA-ESCA.CopyNumberAlteration.Vertex.json.gz \
  -d outputs/tcga/TCGA-ESCA.CopyNumberAlterationOf.Edge.json.gz \
  -d outputs/tcga/TCGA-GBM.CopyNumberAlteration.Vertex.json.gz \
  -d outputs/tcga/TCGA-GBM.CopyNumberAlterationOf.Edge.json.gz \
  -d outputs/tcga/TCGA-HNSC.CopyNumberAlteration.Vertex.json.gz \
  -d outputs/tcga/TCGA-HNSC.CopyNumberAlterationOf.Edge.json.gz \
  -d outputs/tcga/TCGA-KICH.CopyNumberAlteration.Vertex.json.gz \
  -d outputs/tcga/TCGA-KICH.CopyNumberAlterationOf.Edge.json.gz \
  -d outputs/tcga/TCGA-KIRC.CopyNumberAlteration.Vertex.json.gz \
  -d outputs/tcga/TCGA-KIRC.CopyNumberAlterationOf.Edge.json.gz \
  -d outputs/tcga/TCGA-KIRP.CopyNumberAlteration.Vertex.json.gz \
  -d outputs/tcga/TCGA-KIRP.CopyNumberAlterationOf.Edge.json.gz \
  -d outputs/tcga/TCGA-LAML.CopyNumberAlteration.Vertex.json.gz \
  -d outputs/tcga/TCGA-LAML.CopyNumberAlterationOf.Edge.json.gz \
  -d outputs/tcga/TCGA-LGG.CopyNumberAlteration.Vertex.json.gz \
  -d outputs/tcga/TCGA-LGG.CopyNumberAlterationOf.Edge.json.gz \
  -d outputs/tcga/TCGA-LIHC.CopyNumberAlteration.Vertex.json.gz \
  -d outputs/tcga/TCGA-LIHC.CopyNumberAlterationOf.Edge.json.gz \
  -d outputs/tcga/TCGA-LUAD.CopyNumberAlteration.Vertex.json.gz \
  -d outputs/tcga/TCGA-LUAD.CopyNumberAlterationOf.Edge.json.gz \
  -d outputs/tcga/TCGA-LUSC.CopyNumberAlteration.Vertex.json.gz \
  -d outputs/tcga/TCGA-LUSC.CopyNumberAlterationOf.Edge.json.gz \
  -d outputs/tcga/TCGA-MESO.CopyNumberAlteration.Vertex.json.gz \
  -d outputs/tcga/TCGA-MESO.CopyNumberAlterationOf.Edge.json.gz \
  -d outputs/tcga/TCGA-OV.CopyNumberAlteration.Vertex.json.gz \
  -d outputs/tcga/TCGA-OV.CopyNumberAlterationOf.Edge.json.gz \
  -d outputs/tcga/TCGA-PAAD.CopyNumberAlteration.Vertex.json.gz \
  -d outputs/tcga/TCGA-PAAD.CopyNumberAlterationOf.Edge.json.gz \
  -d outputs/tcga/TCGA-PCPG.CopyNumberAlteration.Vertex.json.gz \
  -d outputs/tcga/TCGA-PCPG.CopyNumberAlterationOf.Edge.json.gz \
  -d outputs/tcga/TCGA-PRAD.CopyNumberAlteration.Vertex.json.gz \
  -d outputs/tcga/TCGA-PRAD.CopyNumberAlterationOf.Edge.json.gz \
  -d outputs/tcga/TCGA-READ.CopyNumberAlteration.Vertex.json.gz \
  -d outputs/tcga/TCGA-READ.CopyNumberAlterationOf.Edge.json.gz \
  -d outputs/tcga/TCGA-SARC.CopyNumberAlteration.Vertex.json.gz \
  -d outputs/tcga/TCGA-SARC.CopyNumberAlterationOf.Edge.json.gz \
  -d outputs/tcga/TCGA-SKCM.CopyNumberAlteration.Vertex.json.gz \
  -d outputs/tcga/TCGA-SKCM.CopyNumberAlterationOf.Edge.json.gz \
  -d outputs/tcga/TCGA-STAD.CopyNumberAlteration.Vertex.json.gz \
  -d outputs/tcga/TCGA-STAD.CopyNumberAlterationOf.Edge.json.gz \
  -d outputs/tcga/TCGA-TGCT.CopyNumberAlteration.Vertex.json.gz \
  -d outputs/tcga/TCGA-TGCT.CopyNumberAlterationOf.Edge.json.gz \
  -d outputs/tcga/TCGA-THCA.CopyNumberAlteration.Vertex.json.gz \
  -d outputs/tcga/TCGA-THCA.CopyNumberAlterationOf.Edge.json.gz \
  -d outputs/tcga/TCGA-THYM.CopyNumberAlteration.Vertex.json.gz \
  -d outputs/tcga/TCGA-THYM.CopyNumberAlterationOf.Edge.json.gz \
  -d outputs/tcga/TCGA-UCEC.CopyNumberAlteration.Vertex.json.gz \
  -d outputs/tcga/TCGA-UCEC.CopyNumberAlterationOf.Edge.json.gz \
  -d outputs/tcga/TCGA-UCS.CopyNumberAlteration.Vertex.json.gz \
  -d outputs/tcga/TCGA-UCS.CopyNumberAlterationOf.Edge.json.gz \
  -d outputs/tcga/TCGA-UVM.CopyNumberAlteration.Vertex.json.gz \
  -d outputs/tcga/TCGA-UVM.CopyNumberAlterationOf.Edge.json.gz \
  "python3 transform/dvc/bmeg_file_manifest.py"
```


Sanity check
-------

```
# list all outputs ( except for pubmed text files), make sure they are in some dvc file
ls -1 outputs/**/*.json.gz | grep -v pubmed18 | parallel 'grep  {} *.dvc || echo "no match {}"' | grep 'no match'
```
