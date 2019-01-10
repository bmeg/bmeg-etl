import yaml

filename = 'outputs.bmeg_manifest.dvc'
with open(filename, 'r') as stream:
    dvc = yaml.safe_load(stream)
    with open('outputs/meta/bmeg_file_manifest.txt', 'w+') as fobj:
        fobj.write('\n'.join([d['path'] for d in dvc['deps']]))
    with open('scripts/bmeg_file_manifest.txt', 'w+') as fobj:
        fobj.write('\n'.join([d['path'] for d in dvc['deps']]))
