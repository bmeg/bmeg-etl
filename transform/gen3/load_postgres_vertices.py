import sys

from transform.gen3.gen3_util import get_class_tablename_from_id
from transform.gen3.cli import default_args
from transform.gen3.schemas import make_id


def transform(project_id, vertex_name, output_stream=sys.stdout):
    """Writes transform and load commands to output_stream."""
    table_name = get_class_tablename_from_id(make_id(vertex_name))
    output_stream.write("python transform/gen3/postgres_vertices.py --project_id {} --vertex_name {} | $psql -c \"copy {}(node_id, acl, _sysan,  _props) from stdin  csv delimiter E'\\x01' quote E'\\x02' ;\"".format(project_id, vertex_name, table_name))


if __name__ == "__main__":

    args = default_args()
    transform(args.project_id, args.vertex_name)
