import yaml

""" recreate the command that made the dvc file """

# dvc pipeline show outputs.bmeg_manifest.dvc | grep source
SOURCE_COMMANDS = """
source.ccle.CCLE_DepMap_18q3_maf_20180718.txt.dvc
source.gene_enricher.non_alt_loci_set.json.dvc
source.g2p.all.dvc
source.mc3.v0.2.8.PUBLIC.maf.gz.dvc
source.myvariants.biothings_current_old_hg19.json.gz.dvc
source.myvariant.info.harvested.json.gz.dvc
source.myvariant.info.metadata.fields.json.dvc
source.myvariant.info.metadata.json.dvc
source.ccle.DepMap-2018q3-celllines.csv.dvc
source.ccle.CCLE_DepMap_18q3_RNAseq_RPKM_20180718.gct.dvc
source.ccle.CCLE_NP24.2009_Drug_data_2015.02.24.csv.dvc
source.ccle.CCLE_tpm.tsv.gz.dvc
source.gdsc.GDSC_AUC.csv.dvc
source.ensembl-protein.homo_sapiens.json.dvc
source.ensembl.Homo_sapiens.GRCh37.87.gff3.gz.dvc
source.goa_human.gaf.gz.dvc
source.go.HUMAN_9606_idmapping.dat.gz.dvc
source.go.obo.dvc
source.gtex.GTEx_v7_Annotations_SampleAttributesDS.txt.dvc
source.gtex.GTEx_v7_Annotations_SubjectPhenotypesDS.txt.dvc
source.gtex.GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct.dvc
source.pfam.clans.tsv.dvc
source.pfam.tar.gz.dvc
source.pfam.id_list.txt.dvc
source.tcga.TCGA_ID_MAP.csv.dvc
source.tcga.tcga-genomics.zip.dvc
source.tcga.TCGA_expression_tpm.tsv.gz.dvc
""".strip().split()
# dvc pipeline show outputs.bmeg_manifest.dvc | grep outputs
OUTPUT_COMMANDS = """
outputs.ccle.maf.dvc
outputs.g2p.dvc
outputs.gdc.cases.dvc
outputs.mc3.dvc
outputs.allele.dvc
outputs.ccle.dvc
outputs.ccle.expression.dvc
outputs.ccle.drug_response.dvc
outputs.ccle.expression_tatlow.dvc
outputs.gdsc.dvc
outputs.compound.normalized.dvc
outputs.ensembl-protein.dvc
outputs.ensembl.dvc
outputs.ensembl.missing_transcripts.dvc
outputs.gdc.projects.dvc
outputs.go.gaf2schema.dvc
outputs.go.dvc
outputs.gtex.dvc
outputs.gtex.expression.dvc
outputs.pfam.dvc
outputs.phenotype.dvc
outputs.publication.dvc
outputs.tcga.expression.dvc
outputs.bmeg_manifest.dvc
""".strip().split()


def transform():
    # reverse sort because convention 'source' comes before 'outputs'
    for filename in SOURCE_COMMANDS + OUTPUT_COMMANDS:
        with open(filename, 'r') as stream:
            dvc = yaml.load(stream)
            # not a command, an import
            if 'cmd' not in dvc:
                print('#')
                print('dvc import --file {} --yes \\'.format(filename))
                print(' {}\\'.format(dvc['deps'][0]['path']))
                print(' {}'.format(dvc['outs'][0]['path']))
                continue
            # a dvc run
            print('#')
            print('dvc run --file {} --yes \\'.format(filename))
            for dep in dvc.get('deps', []):
                print('  -o {} \\'.format(dep['path']))
            for out in dvc.get('outs', []):
                print('  -o {} \\'.format(out['path']))
            print('  "{}"'.format(dvc['cmd']))


if __name__ == '__main__':
    transform()
