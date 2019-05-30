import os
import shutil
from transform.gen3.schemas import transform
from yaml import load, FullLoader


def test_schemas():
    """Validates gen3 schema files."""
    this_dir = os.path.dirname(os.path.realpath(__file__))
    outputs_dir = os.path.join(this_dir, 'outputs')
    try:
        shutil.rmtree(outputs_dir)
        os.mkdir(outputs_dir)
    except Exception:
        os.mkdir(outputs_dir)

    schema_files = transform(file_path='source/gen3/edges.json', output_dir=outputs_dir)
    for schema_file in schema_files:
        schema = load(open(schema_file), Loader=FullLoader)
        # cross check
        rp = set(schema['required'])
        sp = set(schema['systemProperties'])
        p = set(schema['properties'].keys())
        link_names = set([l['name'] for l in schema['links']])
        assert sp.issubset(p), 'schema has system property(s) not found in properties {}'.format(sp - p)
        assert rp.issubset(p), 'schema has required property(s) not found in properties {}'.format(rp - p)
        assert link_names.issubset(p), 'schema has link names(s) not found in properties {}'.format(link_names - p)
        # # TODO - once we can load schema and _definitions.yaml
        # from jsonschema.validators import validator_for
        # validator = validator_for(schema)
        # validator.check_schema(schema)
