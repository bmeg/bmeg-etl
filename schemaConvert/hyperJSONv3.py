import click
import os
import pathlib
import glob
import yaml
import json
import string

@click.group()
def cli():
    pass

@cli.command('yaml_dir')
@click.option('--input_path', required=True,
              default='generate-json-schema/coherent.json',
              show_default=True,
              help='Path to schema files directory'
              )
@click.option('--output_path', required=True,
              default='.',
              show_default=True,
              help='Path to cytoscape files directory'
              )
@click.option('--include_links/--exclude_links', default=True)
def yaml_dir(input_path, output_path, include_links):
    """
    Extract a SIF file for import into cytoscape.
    python schema.py cytoscape --input_path aced-bmeg.json --output_path output
    JS expects the resultant schema file to be placed in public/data
    """
    input_path = pathlib.Path(input_path)
    output_path = pathlib.Path(output_path)
    edgeTitles = ["SomaticVariant"]
    copyFields = ["description", "required", "title", "type"]
    assert input_path.is_dir()
    assert output_path.is_dir()
    paths = [file for file in glob.glob(os.path.join(input_path, "*.yaml"))]
    for path in paths:
      if str(path).split("/")[-1][0] != "_":
        with open(os.path.join(output_path, str(path).split("/")[-1]), 'w') as out_file:
            with open(path, "r") as file:
                schema = yaml.safe_load(file)
                new_schema = {}
                new_schema["$schema"] = schema["$schema"]
                new_schema["$id"] = schema["id"]
                for k in schema:
                   if k in copyFields:
                       new_schema[k] = schema[k]
                is_edge = schema["title"] in edgeTitles
                new_schema["links"] = []
                if include_links:
                  new_schema["properties"] = {
                      "links": {
                          "type":"array",
                          "items": {
                              "$ref": "https://json-schema.org/draft/2020-12/links"
                          }
                      }
                  }
                else:
                  new_schema['properties'] = {}
                if "links" in schema:
                  for l in schema["links"]:
                    new_link = {
                        "rel":l["label"],
                        "href":l["target_type"]+"/{id}",
                        "templateRequired": ["id"],
                        "targetSchema": {"$ref": l["target_type"]+'.yaml'},
                        "templatePointers": {'id': "/%s/-/id" % l['label']},
                        "targetHints": {
                            "directionality": ["outbound"],
                            "multiplicity": ["has_one"] if l["multiplicity"].endswith('one') else ["has_many"],
                       }
                    }
                    if is_edge:
                        new_link["targetHints"]["association"] = True
                    else:
                        new_link["targetHints"]["backref"] = schema['id']+('s' if l['multiplicity'].startswith('many') else '')
                    new_schema["links"].append(new_link)
                for p in schema["properties"]:
                    if "targets" in schema["properties"][p] and 'backref' in schema['properties'][p]['targets'][0] and p not in [lin['rel'] for lin in new_schema['links']]:
                        hname = "".join(p.split("_"))[0: -1 if '$ref' in schema["properties"][p] and schema["properties"][p]['$ref'].endswith('many') else None]
                        new_link = {
                            "rel":p,
                            "href":hname+"/{id}",
                            "templateRequired": ["id"],
                            "targetSchema": {'$ref': schema['properties'][p]['targets'][0]['type']['$ref'] if 'targets' in schema['properties'][p] else p[0: -1 if '$ref' in schema["properties"][p] and schema["properties"][p]['$ref'].endswith('many') else None]+'.yaml'},
                            "templatePointers": {'id': "/%s/-/id" % p},
                            "targetHints": {
                                "directionality": ["outbound"],
                                "multiplicity": ["has_many"] if '$ref' in schema["properties"][p] and schema["properties"][p]['$ref'].endswith('many') else ["has_one"],
                            }
                        }
                        if is_edge:
                            new_link["targetHints"]["association"] = True
                        else:
                            try:
                                new_link["targetHints"]["backref"] = schema["properties"][p]["targets"][0]["backref"]
                            except:
                                new_link["targetHints"]["backref"] = schema["id"] + 's'
                        new_schema['links'].append(new_link)
                    elif not('$ref' in schema["properties"][p] and schema['properties'][p]['$ref'].split("/")[-1][0:3] == 'to_'):
                        new_schema["properties"][p] = schema["properties"][p]
                for p in [lin['rel'] for lin in new_schema['links']]:
                    new_schema['properties'][p] = {'type': ['array'], 'items': {'$ref': 'reference.yaml'}}
                out_file.write(yaml.dump(new_schema, sort_keys=False))


if __name__ == '__main__':
  cli()


