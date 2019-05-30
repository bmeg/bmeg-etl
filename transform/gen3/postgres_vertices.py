import uuid
import sys
import os
from bmeg.ioutils import reader
import subprocess
import json
from transform.gen3.cli import default_args


def vertex_files(find_cmd):
    """Gets all Vertex files."""
    p = subprocess.Popen(find_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    p.wait()
    for line in p.stdout.readlines():
        line = line.rstrip().decode('UTF-8')
        if not os.path.isfile(line):
            print('# WARNING: ', line, 'is not a file', file=sys.stderr)
            continue
        yield line


def find_cmd(vertex_name):
    """Formats find command."""
    return 'cat scripts/bmeg_file_manifest.txt | grep Vertex | grep {}'.format(vertex_name)


def transform(project_id, vertex_name, vertex_paths=None, output_stream=sys.stdout, callback=None):
    """Read bmeg json and writes postgres TSV with embedded gen3 json."""
    if not vertex_paths:
        vertex_paths = vertex_files(find_cmd(vertex_name))
    for p in vertex_paths:
        for line in reader(p):
            line = json.loads(line)
            submitter_id = line['gid'].lower()
            node_id = uuid.uuid5(uuid.NAMESPACE_DNS, submitter_id)
            line['data']['project_id'] = project_id
            line['data']['submitter_id'] = submitter_id
            line['node_id'] = node_id
            if callback:
                line = callback(line)
            # copy node_gene(node_id, acl, _sysan,  _props) from stdin  csv delimiter E'\x01' quote E'\x02' ;"
            output_stream.write('{}\x01{}\x01{}\x01{}\n'.format(node_id, '{}', '{}', json.dumps(line['data'], separators=(',', ':'))))

# vertex_name = os.path.basename(p).replace('.Vertex.json', '').replace('.gz','').split('.')[-1]


if __name__ == "__main__":

    args = default_args()
    transform(args.project_id, args.vertex_name)
