import click
import yaml
import glob
import pathlib
import os
from collections import OrderedDict
import yaml.emitter
import yaml.serializer
import yaml.representer
import yaml.resolver


@click.group()
def cli():
    pass


@cli.command('yaml_dir')
@click.option('--input_path', required=True,
              default='schema',
              show_default=True,
              help='Path to schema files directory'
              )
@click.option('--output_path', required=True,
              default='schemaConvert/output',
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
                if "description" not in schema:
                    # aliquot, allele_effect and allele are missing descriptions
                    # was able to get aliquot description from the GDC data dictionary but not sure about the others
                    if schema["$id"] == "Aliquot":
                        schema["description"] = "Pertaining to a portion of the whole. Any one of two or more samples of something, of the same volume or weight"
                    if schema["$id"] == "AlleleEffect":
                        schema["description"] = "TODO: This needs to be filled in with relevent information"
                    if schema["$id"] == "Allele":
                        schema["description"] = "TODO: This needs to be filled in with relevent information"
                elif schema["description"].endswith("\n"):
                    schema["description"] = schema["description"][:-1]

                schema["$id"] = convert_file_to_title(schema["$id"], False)

                if "links" in schema:
                    for i, _ in enumerate(schema["links"]):
                        if "href" in schema["links"][i]:
                            split_values = schema["links"][i]["href"].split("/")
                            schema["links"][i]["href"] = convert_file_to_title(split_values[0], False) + "/" + split_values[1]
                        if "targetSchema" in schema["links"][i] and "$ref" in schema["links"][i]["targetSchema"]:
                            schema["links"][i]["targetSchema"]["$ref"] = convert_file_to_title(schema["links"][i]["targetSchema"]["$ref"], True)
                            if "items" in schema["properties"][schema["links"][i]["rel"]] and "$ref" in schema["properties"][schema["links"][i]["rel"]]["items"]:
                                item = schema["properties"][schema["links"][i]["rel"]]["items"]["$ref"]
                                if item == "reference.yaml" and "eference" not in schema["links"][i]["targetSchema"]["$ref"]:
                                    schema["properties"][schema["links"][i]["rel"]]["items"]["$ref"] = schema["links"][i]["targetSchema"]["$ref"]

                        if "targetHints" in schema["links"][i] and \
                            "backref" in schema["links"][i]["targetHints"] and\
                                isinstance(schema["links"][i]["targetHints"]["backref"], str):
                            schema["links"][i]["targetHints"]["backref"] = [schema["links"][i]["targetHints"]["backref"]]

                if "properties" in schema:
                    for key in schema["properties"].keys():
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
                                del schema["properties"][key]["type"]

                        # Special case row looked like {'oneOf': [{'type': 'null'}, {'$ref': '_definitions.yaml#/genome'}]}
                        elif "oneOf" in schema["properties"][key]:
                            schema["properties"][key]["oneOf"] = schema["properties"][key]["oneOf"][1]
                            # print(schema["properties"][key])
                            # print(schema["properties"][key]["type"], "WHOLE STRUCT: ", schema["properties"][key])

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
