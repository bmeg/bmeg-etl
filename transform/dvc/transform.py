import yaml
import glob
from bmeg.vertex import File, Command
from bmeg.edge import Reads, Writes
from bmeg.util.logging import default_logging
from bmeg.util.cli import default_argument_parser
from bmeg.emitter import new_emitter
""" simple transform of dvc into md """


DEFAULT_DIRECTORY = 'meta'


def transform(
    dvc_path="./*.dvc",
    emitter_name="json",
    output_dir="outputs",
    emitter_directory=DEFAULT_DIRECTORY,
):
    emitter = new_emitter(name=emitter_name, directory=emitter_directory)
    file_keys = ['path', 'md5']
    dups = []
    for filename in glob.iglob(dvc_path, recursive=False):
        with open(filename, 'r') as stream:
            dvc = yaml.load(stream)
            if 'cmd' not in dvc:
                dvc['cmd'] = '# unknown'
            dvc['filename'] = filename
            command = Command(cmd=dvc['cmd'], md5=dvc['md5'], filename=dvc['filename'])
            emitter.emit_vertex(command)
            if 'deps' in dvc:
                for dep in dvc['deps']:
                    if 'md5' not in dep:
                        dep['md5'] = 'stub-md5.{}'.format(dep['path'])
                    for k in [k for k in dep.keys() if k not in file_keys]:
                        del dep[k]
                    file = File(**dep)
                    if file.gid() not in dups:
                        emitter.emit_vertex(file)
                        dups.append(file.gid())
                    emitter.emit_edge(Reads(), command.gid(), file.gid())

            if 'outs' in dvc:
                for out in dvc['outs']:
                    if 'md5' not in out:
                        out['md5'] = 'stub-md5.{}'.format(out['path'])
                    for k in [k for k in out.keys() if k not in file_keys]:
                        del out[k]
                    file = File(**out)
                    if file.gid() not in dups:
                        emitter.emit_vertex(file)
                        dups.append(file.gid())
                    emitter.emit_edge(Writes(), command.gid(), file.gid())
    emitter.close()


if __name__ == '__main__':
    parser = default_argument_parser(prefix_default='g2p')
    parser.add_argument('--dvc_path', type=str,
                        default='./*.dvc',
                        help='Path to dvc files [./*.dvc]')
    # We don't need the first argument, which is the program name
    options = parser.parse_args()
    default_logging(options.loglevel)
    transform(dvc_path=options.dvc_path)
