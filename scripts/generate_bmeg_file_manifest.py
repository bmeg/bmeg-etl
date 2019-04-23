import glob
import os
import subprocess
import shlex
import yaml


files = glob.glob("outputs.*.dvc")

EXCEPTIONS = [
    "outputs/ccle/maf.Allele.Vertex.json.gz",
    "outputs/ccle/maf.Deadletter.Vertex.json.gz",
    "outputs/ccle/drug_response.Compound.Vertex.json.gz",
    "outputs/ccle/drug_response.ResponseTo.Edge.json.gz",
    "outputs/ctrp/ctrp.Compound.Vertex.json.gz",
    "outputs/ctrp/ctrp.ResponseTo.Edge.json.gz",
    "outputs/g2p/Allele.Vertex.json.gz",
    "outputs/g2p/Compound.Vertex.json.gz",
    "outputs/g2p/Deadletter.Vertex.json.gz",
    "outputs/g2p/HasEnvironment.Edge.json.gz",
    "outputs/gdc/Compound.Vertex.json.gz",
    "outputs/gdc/TreatedWith.Edge.json.gz",
    "outputs/gdsc/gdsc.Compound.Vertex.json.gz",
    "outputs/gdsc/gdsc.ResponseTo.Edge.json.gz",
    "outputs/mc3/Allele.Vertex.json.gz",
    "outputs/mc3/Deadletter.Vertex.json.gz",
    "outputs/meta/Command.Vertex.json.gz",
    "outputs/meta/File.Vertex.json.gz",
    "outputs/meta/Reads.Edge.json.gz",
    "outputs/meta/Writes.Edge.json.gz",
    "outputs/meta/bmeg_file_manifest.txt",
    "outputs/tcga/IlluminaHumanMethylation450.Methylation.Vertex.json.gz",
    "outputs/tcga/IlluminaHumanMethylation450.MethylationOf.Edge.json.gz",
    "outputs/tcga/IlluminaHumanMethylation450.MethylationProbe.Vertex.json.gz",
    "outputs/tcga/IlluminaHumanMethylation450.MethylationProbeFor.Edge.json.gz",    
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

for o in sorted(set(outputs)):
    DVC_CMD += " -d {}".format(o)

DVC_CMD += ' "echo generating file manifest..."'
args = shlex.split(DVC_CMD)
subprocess.call(args)
with open('scripts/bmeg_file_manifest.txt', 'w+') as fobj:
    fobj.write('\n'.join(sorted(set(outputs))))
