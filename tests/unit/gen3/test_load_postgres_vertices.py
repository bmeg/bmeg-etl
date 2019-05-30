from transform.gen3.load_postgres_vertices import transform
from io import StringIO


def test_load_vertices():
    """Validates gen3 postgres loader files."""
    output_stream = StringIO()
    transform(vertex_name='Project', project_id='foo', output_stream=output_stream)
    for line in output_stream.getvalue().split('\n'):
        assert '--vertex_name Project' in line, 'Should have a vertex_name'
        assert '--project_id foo' in line, 'Should have a project_id'
        assert 'node_project' in line, 'Should have a table name'
    output_stream.close()
