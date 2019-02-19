import yaml
import glob
import json
from bmeg.util.logging import default_logging
from bmeg.util.cli import default_argument_parser
""" simple transform of dvc into cytoscape """


DEFAULT_DIRECTORY = 'meta'


def transform(
    dvc_path="./*.dvc",
    emitter_name="json",
    output_dir="outputs",
    emitter_directory=DEFAULT_DIRECTORY,
):
    file_keys = ['path', 'md5']
    graph = []
    dups = []
    for filename in glob.iglob(dvc_path, recursive=False):
        with open(filename, 'r') as stream:
            dvc = yaml.load(stream)
            if 'cmd' not in dvc:
                dvc['cmd'] = '# unknown'
            dvc['filename'] = filename
            graph.append({'data': {'id': dvc['md5'], 'cmd': dvc['cmd'], '_label': 'Command'}})
            if 'deps' in dvc:
                for dep in dvc['deps']:
                    if 'md5' not in dep:
                        dep['md5'] = 'stub-md5.{}'.format(dep['path'])
                    for k in [k for k in dep.keys() if k not in file_keys]:
                        del dep[k]
                    file = {'data': {'id': dep['md5'], 'path': dep['path'], '_label': 'File'}}
                    if dep['md5'] not in dups:
                        graph.append(file)
                        dups.append(dep['md5'])
                    graph.append(
                        {'data': {
                            'target': dvc['md5'],
                            'source': dep['md5'],
                            'path': dep['path'],
                            'cmd': dvc['cmd'],
                            '_label': 'READS'
                        }}
                    )
            if 'outs' in dvc:
                for out in dvc['outs']:
                    if 'md5' not in out:
                        out['md5'] = 'stub-md5.{}'.format(out['path'])
                    for k in [k for k in out.keys() if k not in file_keys]:
                        del out[k]
                    file = {'data': {'id': out['md5'], 'path': out['path'], '_label': 'File'}}
                    if out['md5'] not in dups:
                        graph.append(file)
                        dups.append(out['md5'])
                    graph.append(
                        {'data': {
                            'target': out['md5'],
                            'source': dvc['md5'],
                            'path': out['path'],
                            'cmd': dvc['cmd'],
                            '_label': 'WRITES'
                        }}
                    )
    with open('{}/{}/cytoscape.json'.format(output_dir, emitter_directory), 'w') as f:
        f.write('{}\n'.format(json.dumps(graph)))


if __name__ == '__main__':
    parser = default_argument_parser(prefix_default='g2p')
    parser.add_argument('--dvc_path', type=str,
                        default='./*.dvc',
                        help='Path to dvc files [./*.dvc]')
    # We don't need the first argument, which is the program name
    options = parser.parse_args()
    default_logging(options.loglevel)
    transform(dvc_path=options.dvc_path)
