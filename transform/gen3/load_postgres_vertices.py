import sys
import inflection

from transform.gen3.gen3_util import get_class_tablename_from_id
from transform.gen3.cli import default_args


def transform(project_id, vertex_name, output_stream=sys.stdout):
    table_name = get_class_tablename_from_id(inflection.underscore(vertex_name))
    output_stream.write("python transform/gen3/postgres_vertices.py --project_id {} --vertex_name {} | $psql -c \"copy {}(node_id, acl, _sysan,  _props) from stdin  csv delimiter E'\\x01' quote E'\\x02' ;\"".format(project_id, vertex_name, table_name))


if __name__ == "__main__":

    args = default_args()
    transform(args.project_id, args.vertex_name)
