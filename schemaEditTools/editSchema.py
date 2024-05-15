import os
import click
import yaml
import shutil

def copy_files(src_dir, dest_dir):
    try:
        for filename in os.listdir(src_dir):
            src_file = os.path.join(src_dir, filename)
            dest_file = os.path.join(dest_dir, filename)
            shutil.copy(src_file, dest_file)
    except OSError as e:
        print(f"Error: {e.strerror}")

def check_and_remove_dir(dest):
    if os.path.isdir(dest):
        for file in os.listdir(dest):
            rem_file = os.path.join(dest,file)
            os.unlink(rem_file)
        os.rmdir(dest)

def copy_properties(schema_dir, source_schema_name, copy_schema_name):
    with open(f"{os.path.join(schema_dir, copy_schema_name)}.yaml", 'r') as copy_file:
        copy_schema_data = yaml.safe_load(copy_file)

    with open(f"{os.path.join(schema_dir, source_schema_name)}.yaml", 'r') as source_file:
        source_schema_data = yaml.safe_load(source_file)

    source_properties = source_schema_data["properties"]
    for key, value in source_properties.items():
        if key in copy_schema_data["properties"] and key not in ["id", "resourceType"]:
            print(f"key {key} already exists in dest {copy_schema_data["properties"].keys()}")
            raise Exception
        if key not in ["id", "resourceType"]:
            copy_schema_data["properties"][key] = value

    output_directory = schema_dir + "_name"
    os.makedirs(output_directory, exist_ok=True)
    with open(os.path.join(output_directory, f"{copy_schema_name}.yaml"), 'w') as file:
        yaml.dump(copy_schema_data, file, default_flow_style=False)


def remove_link(input_dir, output_dir, schema_name, link_name):
    with open(f"{os.path.join(input_dir, schema_name)}.yaml", 'r') as file:
        schema_data = yaml.safe_load(file)

    schema_data['links'] = [link for link in schema_data['links'] if link.get('rel').split("_")[0] != link_name]

    if 'properties' in schema_data:
        schema_data['properties'].pop(link_name, None)

    os.makedirs(output_dir, exist_ok=True)

    with open(os.path.join(output_dir, f"{schema_name}.yaml"), 'w') as file:
        yaml.dump(schema_data, file, default_flow_style=False)

    click.echo(f"Link '{link_name}' and property '{link_name}' removed. Modified schema saved to '{output_dir}/'.")

def add_link(input_dir, output_dir, schema_name, link_name, href, directionality, multiplicity, type):
    with open(f"{os.path.join(input_dir, schema_name)}.yaml", 'r') as file:
        schema_data = yaml.safe_load(file)
    backref =  f"{link_name}_{schema_data["$id"].lower()}"
    schema_data['links'].append({
        "rel":f"{link_name}_{href}",
        "href": f"{href}/" + '{' + 'id' + '}',
        "templateRequired": ["id"],
        "targetSchema": {"$ref": f"{href}.yaml"},
        "templatePointers": {"id":f"/{link_name}/-/id"},
        "targetHints": {"directionality": [directionality],
                        "multiplicity": [multiplicity],
                        "backref": [backref]
                       }
    })

    schema_data['properties'][link_name] = {
        "element_property":True,
        "type":type,
        "backref":backref,
        "items": { "$ref": f"{href}.yaml"}
    }


    os.makedirs(output_dir, exist_ok=True)

    with open(os.path.join(output_dir, f"{schema_name}.yaml"), 'w') as file:
        yaml.dump(schema_data, file, default_flow_style=False)

    click.echo(f"Link '{link_name}' and property '{link_name}' added. Modified schema saved to '{output_dir}/'.")

if __name__ == '__main__':
    # remove sample, command, file, case, compound, aliquot, program, project

    iceberg = "schemaConvert/iceberg"
    revised = "schemaConvert/revisedSchemas"
    schema_lifted = "schemaConvert/schema_lifted"
    dest = "schemaConvert/ACED_BMEG_UNIFIED"

    # Add and remove links from old graph
    add_link(iceberg, schema_lifted, "Task", "specimen", "Specimen", "outbound", "has_many", "array")
    add_link(schema_lifted, schema_lifted, "Task", "documents", "DocumentReference", "outbound", "has_many", "array")

    add_link(iceberg, schema_lifted, "Specimen", "phenotypes", "Phenotype", "outbound", "has_many", "array")
    add_link(schema_lifted, schema_lifted, "Specimen", "projects", "ResearchStudy", "outbound", "has_many", "array")
    add_link(schema_lifted, schema_lifted, "Specimen", "specimen", "DocumentReference", "outbound", "has_many", "array")

    add_link(iceberg, schema_lifted, "Substance", "projects", "ResearchStudy", "outbound", "has_many", "array")
    add_link(iceberg, schema_lifted, "Patient", "substances", "Substance", "outbound", "has_many", "array")
    add_link(revised, schema_lifted, "Phenotype", "patients", "Patient",  "outbound", "has_many", "array")

    add_link(revised, schema_lifted, "GenePhenotypeAssociation", "substances", "Substance", "outbound", "has_many", "array")
    remove_link(schema_lifted, schema_lifted, "GenePhenotypeAssociation", "compounds")

    add_link(revised, schema_lifted, "ProteinCompoundAssociation", "substance", "Substance", "outbound", "has_one", "array")
    remove_link(schema_lifted, schema_lifted, "ProteinCompoundAssociation", "compound")

    add_link(revised, schema_lifted, "DrugResponse", "substances", "Substance", "outbound", "has_many", "array")
    remove_link(schema_lifted, schema_lifted, "DrugResponse", "compounds")
    remove_link(schema_lifted, schema_lifted, "DrugResponse", "aliquot")

    add_link(revised, schema_lifted, "CopyNumberAlteration", "specimen", "Specimen", "outbound", "has_one", "array")
    remove_link(schema_lifted, schema_lifted, "CopyNumberAlteration", "aliquot")

    add_link(revised, schema_lifted, "SomaticCallset", "specimen", "Specimen", "outbound", "has_many", "array")
    remove_link(schema_lifted, schema_lifted, "SomaticCallset", "aliquots")

    add_link(revised, schema_lifted, "GeneExpression", "specimen", "Specimen", "outbound", "has_one", "array")
    remove_link(schema_lifted, schema_lifted, "GeneExpression", "aliquot")

    add_link(revised, schema_lifted, "Methylation", "specimen", "Specimen", "outbound", "has_one", "array")
    remove_link(schema_lifted, schema_lifted, "Methylation", "aliquot")

    add_link(revised, schema_lifted, "TranscriptExpression", "specimen", "Specimen", "outbound", "has_one", "array")
    remove_link(schema_lifted, schema_lifted, "TranscriptExpression", "aliquot")

    # Remove nodes from schema dir
    node_names = "Sample Command File Case Compound Aliquot Program Project".split(" ")
    node_files = [s + ".yaml" for s in node_names]
    for file in node_files:
        rem_file = os.path.join("schemaConvert/revisedSchemas",file)
        if os.path.isfile(rem_file):
            os.unlink(rem_file)

    # Check for output directory from a previous run and remove it if it exists
    check_and_remove_dir(dest)
    os.mkdir(dest)

    # Copy files from iceberg, and revised to output dir
    copy_files(revised, dest)
    copy_files(iceberg, dest)

    # Copy and replace files from intermediate dir to output dir
    for file in [file for file in os.listdir(schema_lifted) if str(file).endswith(".yaml")]:
        os.replace(os.path.join(schema_lifted, file), os.path.join(dest, file))

    # Check for staging dir from a previous run and remove if it exists
    check_and_remove_dir(schema_lifted)
    check_and_remove_dir(revised)









