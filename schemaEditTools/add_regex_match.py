import os
import yaml


if __name__ == '__main__':
    schema_dir = '../schema'
    schema_files = [f for f in os.listdir(schema_dir) if str(f).endswith(".yaml")]
    schemaNames = [scf.rstrip(".yaml") for scf in schema_files]
    for schema in schema_files:
        with open(f"{schema_dir}/{schema}", "r")as rf:
            schema_data = yaml.safe_load(rf)
            for link in schema_data["links"]:
                if link["href"] == "Resource/{id}":
                    with open(f"{schema}", "w") as wf:
                        for sch in schemaNames:
                            new_link = {
                                "href": "id",
                                "rel": link["rel"] + f"_{sch}",
                                "targetHints": {
                                    "backref": [
                                        link["rel"] + f"_{schema_data["$id"].lower()}"
                                    ],
                                    "direction":
                                        list(link["targetHints"]["direction"])
                                    ,
                                    "multiplicity": [
                                        'has_one'
                                    ],
                                    "regex_match":[
                                        f"{sch}/*"
                                    ]
                                },
                                "targetSchema": {
                                    "$ref": f"{sch}.yaml"
                                },
                                "templatePointers":{
                                    "id": link["templatePointers"]["id"]
                                },
                                "templateRequired":[
                                    "id"
                                ]
                            }
                            schema_data["links"].append(new_link)
                        yaml.dump(schema_data, wf)
