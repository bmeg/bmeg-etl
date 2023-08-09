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
    links = {}
    schemas = {}

    for path in paths:
        if str(path).split("/")[-1][0] != "_":
            with open(path, "r") as file:
                schema = yaml.safe_load(file)
                new_schema = {}
                new_schema["$schema"] = schema["$schema"]
                new_schema["$id"] = schema["id"]
                links[schema['id']] = {}
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
                #case: link in 'links' field: highest priority
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
                    links[schema['id']][l['target_type']] = {
                        'priority': 4,
                        'link': new_link
                    }
                    #new_schema["links"].append(new_link)
                for p in schema["properties"]:
                    #links in properties get a priority of 1 that is added onto depending on the link properties
                    if "targets" in schema["properties"][p] or ('$ref' in schema["properties"][p] and schema['properties'][p]['$ref'].split("/")[-1][0:3] == 'to_'):
                      priority = 1
                      
                      hname = "".join(p.split("_"))[0: -1 if '$ref' in schema["properties"][p] and schema["properties"][p]['$ref'].endswith('many') else None]
                      new_link = {
                          "rel":p,
                          "href":hname+"/{id}",
                          "templateRequired": ["id"],
                          "targetSchema": {'$ref': schema['properties'][p]['targets'][0]['type']['$ref'] if 'targets' in schema['properties'][p] else p[0: -1 if '$ref' in schema["properties"][p] and schema["properties"][p]['$ref'].endswith('many') else None]+'.yaml'},
                          "templatePointers": {'id': "/%s/-/id" % p},
                          "targetHints": {
                              "directionality": ["outbound"],
                              "multiplicity": ["has_one"],
                          }
                      }

                      if is_edge:
                        #highest priority if the object is an edge object
                        priority += 3
                        new_link['targetHints']['association'] = True
                      else:
                        if 'targets' in schema['properties'][p] and 'backref' in schema['properties'][p]['targets'][0]:
                          new_link['targetHints']['backref'] = schema['properties'][p]['targets'][0]['backref']
                        else:
                          new_link['targetHints']['backref'] = schema['id']+'s'
                      if 'targets' in schema['properties'][p]:
                        priority += 1
                      if '$ref' in schema["properties"][p] and schema['properties'][p]['$ref'].endswith('one'):
                        priority += 1
                      

                      if new_link['targetSchema']['$ref'].split('.')[0] not in links[new_schema["$id"]]:
                        links[new_schema['$id']][new_link['targetSchema']['$ref'].split('.')[0]] = {'priority': priority, 'link': new_link}
                      #  new_schema['properties'][p] = {'type': ['array'], 'items': {'$ref': 'reference.yaml'}}
                      #elif p in [lin['rel'] for lin in new_schema['links']]:
                      #  new_schema['properties'][p] = {'type': ['array'], 'items': {'$ref': 'reference.yaml'}}

                    else:
                        new_schema["properties"][p] = schema["properties"][p]
                schemas[schema['id']] = new_schema
               #out_file.write(yaml.dump(new_schema, sort_keys=False))
    for v1 in links:
      for v2 in links[v1]:
        if v2 in links: 
          if v1 in links[v2]:
            if 'priority' in links[v1][v2]:
              priorityDiff = links[v1][v2]['priority'] - links[v2][v1]['priority']
              if priorityDiff == 0 and v1 != v2:
                print('duplicate links: %s and %s' % (v1, v2))
              else:
                a = v1 if priorityDiff > 0 else v2
                b = v2 if priorityDiff > 0 else v1
                schemas[a]['links'].append(links[a][b]['link'])
                schemas[a]['properties'][links[a][b]['link']['rel']] = {'type': ['array'], 'items': {'$ref': 'reference.yaml'}}
              links[v2][v1].pop('priority')
          else:
            a = v1
            b = v2
            schemas[a]['links'].append(links[a][b]['link'])
            schemas[a]['properties'][links[a][b]['link']['rel']] = {'type': ['array'], 'items': {'$ref': 'reference.yaml'}}
        else:
           print(v1 + ': '+ v2 + ' dropped')

    for path in paths:
      if str(path).split("/")[-1][0] != "_":
        with open(os.path.join(output_path, str(path).split("/")[-1]), 'w') as out_file:
            out_file.write(yaml.dump(schemas[str(path).split("/")[-1].split('.')[0]], sort_keys=False))

if __name__ == '__main__':
  cli()
