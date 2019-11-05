import yaml
import sys

""" recreate the command that made the dvc file """

# pipe the following into std in
# dvc pipeline show outputs.bmeg_manifest.dvc | grep source
# dvc pipeline show outputs.bmeg_manifest.dvc | grep outputs


def transform():
    # reverse sort because convention 'source' comes before 'outputs'
    for filename in sys.stdin:
        filename = filename.strip()
        with open(filename, 'r') as stream:
            dvc = yaml.safe_load(stream)
            # not a command, an import
            if 'cmd' not in dvc:
                print('#')
                print('dvc import --file {} --overwrite-dvcfile -w ../.. \\'.format(filename))
                print(' {}\\'.format(dvc['deps'][0]['path']))
                print(' {}'.format(dvc['outs'][0]['path']))
                continue
            # a dvc run
            print('#')
            print('dvc run --file {} --overwrite-dvcfile -w ../.. \\'.format(filename))
            for dep in dvc.get('deps', []):
                print('  -d {} \\'.format(dep['path']))
            for out in dvc.get('outs', []):
                print('  -o {} \\'.format(out['path']))
            if 'cmd' in dvc:
                print('  "{}"'.format(dvc['cmd']))


if __name__ == '__main__':
    transform()
