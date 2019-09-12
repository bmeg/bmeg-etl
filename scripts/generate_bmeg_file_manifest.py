import glob
import os
import subprocess
import shlex
import yaml


files = glob.glob("outputs/**/*.dvc")

EXCEPTIONS = [
    # unnormalized Compounds
    "outputs/ccle/drug_response.Compound.Vertex.json.gz",
    "outputs/ctrp/drug_response.Compound.Vertex.json.gz",
    "outputs/gdsc/drug_response.Compound.Vertex.json.gz",
    "outputs/g2p/Compound.Vertex.json.gz",
    "outputs/gdc/Compound.Vertex.json.gz",
    "outputs/dgidb/Compound.Vertex.json.gz",
    # unnormalized Compound edges
    "outputs/ccle/drug_response.DrugResponse_Compounds_Compound.Edge.json.gz",
    "outputs/ccle/drug_response.Compound_DrugResponses_DrugResponse.Edge.json.gz",
    "outputs/ctrp/drug_response.DrugResponse_Compounds_Compound.Edge.json.gz",
    "outputs/ctrp/drug_response.Compound_DrugResponses_DrugResponse.Edge.json.gz",
    "outputs/gdsc/drug_response.DrugResponse_Compounds_Compound.Edge.json.gz",
    "outputs/gdsc/drug_response.Compound_DrugResponses_DrugResponse.Edge.json.gz",
    "outputs/g2p/G2PAssociation_Compounds_Compound.Edge.json.gz",
    "outputs/g2p/Compound_G2PAssociations_G2PAssociation.Edge.json.gz",
    "outputs/gdc/Case_Compounds_Compound.Edge.json.gz",
    "outputs/gdc/Compound_Cases_Case.Edge.json.gz",
    "outputs/gdc/Compound_Projects_Project.Edge.json.gz ",
    "outputs/gdc/Project_Compounds_Compound.Edge.json.gz",
    "outputs/ccle/drug_response.Project_Compounds_Compound.Edge.json.gz",
    "outputs/ccle/drug_response.Compound_Projects_Project.Edge.json.gz",
    "outputs/ctrp/drug_response.Project_Compounds_Compound.Edge.json.gz",
    "outputs/ctrp/drug_response.Compound_Projects_Project.Edge.json.gz",
    "outputs/gdsc/drug_response.Project_Compounds_Compound.Edge.json.gz",
    "outputs/gdsc/drug_response.Compound_Projects_Project.Edge.json.gz",
    "outputs/dgidb/G2PAssociation_Compounds_Compound.Edge.json.gz",
    "outputs/dgidb/Compound_G2PAssociations_G2PAssociation.Edge.json.gz",
    # unnormalized Phenotypes
    "outputs/ccle/ccle.Phenotype.Vertex.json.gz",
    "outputs/ctrp/ctrp.Phenotype.Vertex.json.gz",
    "outputs/g2p/Phenotype.Vertex.json.gz",
    "outputs/gdsc/gdsc.Phenotype.Vertex.json.gz",
    # unnormalized Phenotype edges
    "outputs/ccle/ccle.Case_Phenotypes_Phenotype.Edge.json.gz",
    "outputs/ccle/ccle.Phenotype_Cases_Case.Edge.json.gz",
    "outputs/ccle/ccle.Sample_Phenotypes_Phenotype.Edge.json.gz",
    "outputs/ccle/ccle.Phenotype_Samples_Sample.Edge.json.gz",
    "outputs/ctrp/ctrp.Case_Phenotypes_Phenotype.Edge.json.gz",
    "outputs/ctrp/ctrp.Phenotype_Cases_Case.Edge.json.gz",
    "outputs/ctrp/ctrp.Sample_Phenotypes_Phenotype.Edge.json.gz",
    "outputs/ctrp/ctrp.Phenotype_Samples_Sample.Edge.json.gz",
    "outputs/gdsc/gdsc.Case_Phenotypes_Phenotype.Edge.json.gz",
    "outputs/gdsc/gdsc.Phenotype_Cases_Case.Edge.json.gz",
    "outputs/gdsc/gdsc.Sample_Phenotypes_Phenotype.Edge.json.gz",
    "outputs/gdsc/gdsc.Phenotype_Samples_Sample.Edge.json.gz",
    "outputs/g2p/G2PAssociation_Phenotypes_Phenotype.Edge.json.gz",
    "outputs/g2p/Phenotype_G2PAssociations_G2PAssociation.Edge.json.gz",
    "outputs/gdc/Case_Phenotypes_Phenotype.Edge.json.gz",
    "outputs/gdc/Phenotype_Cases_Case.Edge.json.gz",
    "outputs/gdc/Sample_Phenotypes_Phenotype.Edge.json.gz",
    "outputs/gdc/Phenotype_Samples_Sample.Edge.json.gz",
    # Deadletter
    "outputs/g2p/Deadletter.Vertex.json.gz",
    "outputs/mc3/Deadletter.Vertex.json.gz",
    # unnormalized Allele
    "outputs/ccle/maf.Allele.Vertex.json.gz",
    "outputs/g2p/Allele.Vertex.json.gz",
    "outputs/mc3/Allele.Vertex.json.gz",
    # Meta files
    "outputs/meta/Command.Vertex.json.gz",
    "outputs/meta/File.Vertex.json.gz",
    "outputs/meta/Command_Reads_File.json.gz",
    "outputs/meta/File_InputTo_Command.json.gz",
    "outputs/meta/Command_Writes_File.Edge.json.gz",
    "outputs/meta/File_CreatedBy_Command.json.gz",
    "outputs/meta/bmeg_file_manifest.txt",
    # Methylation
    "outputs/tcga/IlluminaHumanMethylation450.Methylation.Vertex.json.gz",
    "outputs/tcga/IlluminaHumanMethylation450.MethylationProbe.Vertex.json.gz",
    "outputs/tcga/IlluminaHumanMethylation450.Aliquot_Methylations_Methylation.Edge.json.gz",
    "outputs/tcga/IlluminaHumanMethylation450.Methylation_Aliquot_Aliquot.Edge.json.gz",
    "outputs/tcga/IlluminaHumanMethylation450.MethylationProbe_Gene_Gene.Edge.json.gz",
    "outputs/tcga/IlluminaHumanMethylation450.Gene_MethylationProbes_MethylationProbe.Edge.json.gz"
]

print("generating DVC command...")

DVC_CMD = "dvc run --file outputs.bmeg_manifest.dvc --yes --ignore-build-cache"

outputs = []
for f in files:
    with open(f, "r") as stream:
        dvc = yaml.safe_load(stream)
        if "outs" not in dvc:
            print(f, "has no outputs...")
            continue
        for d in dvc["outs"]:
            if d["path"] in EXCEPTIONS:
                print("excluding {}...".format(d["path"]))
                continue
            if os.path.isfile(d["path"]):
                outputs.append(d["path"])
            elif os.path.isdir(d["path"]):
                ofiles = glob.glob(os.path.join(d["path"], "**", "*.Vertex.json.gz"), recursive=True) + glob.glob(os.path.join(d["path"], "**", "*.Edge.json.gz"), recursive=True)
                for of in ofiles:
                    if of in EXCEPTIONS:
                        print("excluding {}...".format(of))
                        continue
                    outputs.append(of)

final_outputs = []
for o in sorted(set(outputs)):
    if not (o.endswith(".Edge.json.gz") or o.endswith(".Vertex.json.gz")):
        print("excluding {}...".format(o))
        continue
    DVC_CMD += " -d {}".format(o)
    final_outputs.append(o)

DVC_CMD += ' "echo generating file manifest..."'
args = shlex.split(DVC_CMD)
subprocess.call(args)
with open('scripts/bmeg_file_manifest.txt', 'w+') as fobj:
    fobj.write('\n'.join(final_outputs))
