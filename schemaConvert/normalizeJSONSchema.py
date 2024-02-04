import click
import yaml
import glob
import pathlib
import os
import re
from collections import OrderedDict
import yaml.emitter
import yaml.serializer
import yaml.representer
import yaml.resolver


@click.group()
def cli():
    pass


def camelCase(st):
    output = ''.join(x for x in st.title() if x.isalnum())
    out = output[1:]
    # fix special case
    if "Datetime" in str(output[1:]):
        out = str(output[1:]).strip("Datetime") + "DateTime"
    return output[0].lower() + out

@cli.command('yaml_dir')
@click.option('--input_path', required=True,
              default='schema',
              show_default=True,
              help='Path to schema files directory'
              )
@click.option('--output_path', required=True,
              default='schemaConvert/revisedSchemas',
              show_default=True,
              help='Path to cytoscape files directory'
              )
def yaml_dir(input_path, output_path):
    input_path = pathlib.Path(input_path)
    output_path = pathlib.Path(output_path)
    assert input_path.is_dir()
    assert output_path.is_dir()
    paths = [file for file in glob.glob(os.path.join(input_path, "*.yaml"))]
    schemas = {}
    for path in paths:
        if str(path).split("/")[-1][0] != "_":
            with open(path, "r") as file:
                schema = yaml.safe_load(file)
                if "required" in schema:
                    schema["required"] = [camelCase(elem) for elem in schema["required"] if "_" in elem]
                if "properties" in schema:
                    for key in list(schema["properties"].keys()):
                        schema["properties"][camelCase(key)] = schema["properties"][key]
                        if key != camelCase(key):
                            del schema["properties"][key]
                        key = camelCase(key)

                        if "type" in schema["properties"][key]:
                            prop_type_field = schema["properties"][key]["type"]
                            # Logic to simplify the property typing to match FHIR conventions
                            if "type" in schema["properties"][key] and isinstance(prop_type_field, list):
                                if len(prop_type_field) == 1 and prop_type_field[0] not in ["null", "object"]:
                                    schema["properties"][key]["type"] = schema["properties"][key]["type"][0]
                                elif len(prop_type_field) == 2 and prop_type_field[0] in ["null", "object"] and prop_type_field[1] not in ["null", "object"]:
                                    schema["properties"][key]["type"] = schema["properties"][key]["type"][1]
                                elif len(prop_type_field) == 2 and prop_type_field[1] in ["null", "object"] and prop_type_field[0] not in ["null", "object"]:
                                    schema["properties"][key]["type"] = schema["properties"][key]["type"][0]
                                else:
                                    del schema["properties"][key]["type"]

                        # Special case row looked like {'oneOf': [{'type': 'null'}, {'$ref': '_definitions.yaml#/genome'}]}
                        elif "oneOf" in schema["properties"][key]:
                            schema["properties"][key]["oneOf"][0]["type"] = "array"

                        # Everything in BMEG is an element property. Not everything in FHIR is.  ex: fhir_comments propert in DocumentRference
                        schema["properties"][key]["element_property"] = True

                    # resourceType property is very useful in ETL knowing what the node type is
                    # without having to rely on the file name. Should probably be added.
                    schema["properties"]["resourceType"] = {"default": schema["$id"],
                                                            "type": "string",
                                                            "description": "One of the resource types defined as part of BMEG"
                                                            }
                if "description" not in schema:
                    # aliquot, allele_effect and allele are missing descriptions
                    # was able to get aliquot description from the GDC data dictionary but not sure about the others
                    if schema["$id"] == "Aliquot":
                        schema["description"] = "Pertaining to a portion of the whole. Any one of two or more samples of something, of the same volume or weight."
                    if schema["$id"] == "AlleleEffect":
                        schema["description"] = "TODO: This needs to be filled in with relevent information"
                    if schema["$id"] == "Allele":
                        schema["description"] = "TODO: This needs to be filled in with relevent information"
                elif schema["description"].endswith("\n"):
                    schema["description"] = schema["description"][:-1]

                schema["$id"] = convert_file_to_title(schema["$id"], False)

                if "links" in schema:
                    for i, _ in enumerate(schema["links"]):
                        schema_links_i = schema["links"][i]
                        if "href" in schema_links_i:
                            split_values = schema_links_i["href"].split("/")
                            schema_links_i["href"] = convert_file_to_title(split_values[0], False) + "/" + split_values[1]
                        if "targetSchema" in schema_links_i and "$ref" in schema_links_i["targetSchema"]:
                            schema_links_i_rel_camel_case = camelCase(schema_links_i["rel"])
                            schema_links_i["targetSchema"]["$ref"] = convert_file_to_title(schema_links_i["targetSchema"]["$ref"], True)
                            if "items" in schema["properties"][schema_links_i_rel_camel_case] and "$ref" in schema["properties"][schema_links_i_rel_camel_case]["items"]:
                                item = schema["properties"][schema_links_i_rel_camel_case]["items"]["$ref"]
                                if item == "reference.yaml" and "eference" not in schema_links_i["targetSchema"]["$ref"]:
                                    schema["properties"][schema_links_i_rel_camel_case]["items"]["$ref"] = schema_links_i["targetSchema"]["$ref"]

                        if "templatePointers" in schema_links_i:
                            string_ex = schema_links_i["templatePointers"]["id"].split("/")
                            camel_case = camelCase(string_ex[1])
                            if "_" in string_ex[1]:
                                # This implies that only "/projects/-/id" exist and not /projects/-/name/-/id
                                schema_links_i["templatePointers"]["id"] = "/" + camel_case + "/" + "/".join(string_ex[2:])
                            if "rel" in schema_links_i:
                                schema_links_i["rel"] = camel_case + "_" + schema_links_i["targetSchema"]["$ref"].strip(".yaml")

                        if "targetHints" in schema_links_i and \
                            "backref" in schema_links_i["targetHints"] and\
                                isinstance(schema_links_i["targetHints"]["backref"], str):
                            schema["links"][i]["targetHints"]["backref"] = [schema_links_i["targetHints"]["backref"]]
                        if "backref" in schema_links_i['targetHints']:
                            for i, v in enumerate(schema_links_i['targetHints']["backref"]):
                                word = v
                                # Get rid of pluralization of backref since that style doesn't exist in fhir
                                if word.endswith("s"):
                                    word = word[:-1]
                                # Another one that assumes that there aren't any /projects/-/name/-/id refs in Bmeg
                                schema_links_i['targetHints']["backref"][i] = schema_links_i["rel"].split("_")[0] + "_"+ word

                        if "targetHints" in schema_links_i and "backref" in schema_links_i["targetHints"]:
                            schema["properties"][schema_links_i_rel_camel_case]["backref"] = schema_links_i["targetHints"]["backref"][0]

        # Reorder schema keys so that new keys that are added aren't appended to the end of the file
        order = ["$schema", "$id", "title", "type", "description", "required", "links", "properties"]
        ordered_schema = dict(OrderedDict(sorted(schema.items(), key=lambda x: order.index(x[0]))))
        schemas[schema['$id']] = ordered_schema

    for path in paths:
        if str(path).split("/")[-1][0] != "_":
            file_name = convert_file_to_title(str(path).split("/")[-1], True)
            with open(os.path.join(output_path, file_name), 'w') as out_file:
                out_file.write(yaml.dump(schemas[file_name.split(".yaml")[0]], Dumper=PrettyDumper))


"""
    YAML pretty printer to match the YAML style of the FHIR yaml json schema
    Taken from https://github.com/nexB/saneyaml/blob/main/src/saneyaml.py
    and https://stackoverflow.com/questions/60226848/array-does-not-have-indent-or-space-in-pyyaml
"""
class IndentingEmitter(yaml.emitter.Emitter):
    def increase_indent(self, flow=False, indentless=False,):
        """Ensure that lists items are always indented."""
        return super().increase_indent(
            flow=False,
            indentless=False,
        )


class PrettyDumper(
    IndentingEmitter,
    yaml.serializer.Serializer,
    yaml.representer.Representer,
    yaml.resolver.Resolver,
):
    def __init__(
        self,
        stream,
        default_style=None,
        default_flow_style=False,
        canonical=None,
        indent=None,
        width=None,
        allow_unicode=None,
        line_break=None,
        encoding=None,
        explicit_start=None,
        explicit_end=None,
        version=None,
        tags=None,
        sort_keys=False
    ):
        IndentingEmitter.__init__(
            self,
            stream,
            canonical=canonical,
            indent=indent,
            width=width,
            allow_unicode=allow_unicode,
            line_break=line_break,
        )
        yaml.serializer.Serializer.__init__(
            self,
            encoding=encoding,
            explicit_start=explicit_start,
            explicit_end=explicit_end,
            version=version,
            tags=tags,
        )
        yaml.representer.Representer.__init__(
            self,
            default_style=default_style,
            default_flow_style=default_flow_style,
            sort_keys=False,
        )
        yaml.resolver.Resolver.__init__(self)


def convert_file_to_title(file, path):
    file = file.split('.yaml')[0]
    if "_" in file:
        file = file.split("_")
        for i, _ in enumerate(file):
            file[i] = file[i].capitalize()
        file = "".join(file)
    elif not file[0].isupper():
        file = file.capitalize()

    if path:
        file = file + ".yaml"
    return file

if __name__ == '__main__':
    cli()
