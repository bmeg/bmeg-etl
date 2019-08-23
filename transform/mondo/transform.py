import json
import pydash

from bmeg import Phenotype, Project, Phenotype_ParentTerms_Phenotype
from bmeg.emitter import JSONEmitter


def transform(input_file="source/mondo/mondo.json",
              emitter_prefix=None,
              emitter_directory="mondo"):

    emitter = JSONEmitter(directory=emitter_directory, prefix=emitter_prefix)

    with open(input_file, "r") as json_file:
        data = json.load(json_file)
        for term in data["graphs"][0]["nodes"]:
            if term.get("type", "UNKNOWN") != "CLASS":
                continue
            if pydash.get(term, "meta.deprecated", False):
                continue
            # should we include an array of exactMatch synonyms is available?
            pheno = Phenotype(
                term=term.get("lbl", ""),
                term_id=term["id"].replace("http://purl.obolibrary.org/obo/", ""),
                id=Phenotype.make_gid(term["id"].replace("http://purl.obolibrary.org/obo/", "")),
                project_id=Project.make_gid("Reference")
            )
            emitter.emit_vertex(pheno)
        for rel in data["graphs"][0]["edges"]:
            if rel["pred"] != "is_a":
                continue
            emitter.emit_edge(
                Phenotype_ParentTerms_Phenotype(
                    from_gid=Phenotype.make_gid(rel["sub"].replace("http://purl.obolibrary.org/obo/", "")),
                    to_gid=Phenotype.make_gid(rel["obj"].replace("http://purl.obolibrary.org/obo/", ""))
                ),
                emit_backref=True
            )

    emitter.close()


if __name__ == "__main__":
    transform()
