import os
from transform.gen3.postgres_vertices import transform
from yaml import load, FullLoader
from io import StringIO
from jsonschema import validate
import json


def test_load_vertices():
    """Validates gen3 postgres loader files."""
    this_dir = os.path.dirname(os.path.realpath(__file__))

    schema_file = os.path.join(this_dir, 'project.yaml')
    project = load(open(schema_file), Loader=FullLoader)

    output_stream = StringIO()
    transform('ohsu-test', 'Project', vertex_paths=None, output_stream=output_stream)
    for line in output_stream.getvalue().split('\n'):
        if len(line) < 1:
            continue
        node_id, acl, _sysan, _props = line.split('\x01')
        assert node_id, 'Should have a node_id'
        assert acl, 'Should have a acl'
        assert _sysan, 'Should have a _sysan'
        assert _props, 'Should have a _props'
        _props = json.loads(_props)
        # {"project_id":"ohsu-test","gdc_attributes":{},"submitter_id":"project:gtex"}
        validate(_props, project)
    output_stream.close()
