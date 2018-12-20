import yaml
import glob
import types
import bmeg.ioutils
from tabulate import tabulate

""" simple transform of dvc into md """


DEFAULT_DIRECTORY = 'meta'

HUGO_HEADER = """
---
title: Data Sources
sidebar: true
menu:
  main:
    parent: Building Graph
    weight: 2
---
"""

def cardify(papers):
    """create a 'card' for each paper."""
    papers.sort(key=lambda item: (len(item['comment']), item))
    lines = []
    lines.append('{{% card-container %}}')
    for paper in papers:
        source = paper['path'].replace('source/','').replace('/*','')
        lines.append('{{{{% card title="{}" %}}}}'.format(source))
        lines.append(paper['comment'])
        for publication in paper['publications'].split(','):
            lines.append('[paper]({})'.format(publication))
        lines.append('{{% /card %}}')

    lines.append('{{% /card-container %}}')
    return '\n'.join(lines)


def transform(
    source_dir=".",
    papers_path='source/meta/bmeg-source-data-papers.tsv',
    emitter_directory=None
):
    papers_reader = bmeg.ioutils.read_tsv(papers_path)
    papers = []
    for source in papers_reader:
        papers.append(source)

    if not emitter_directory:
        emitter_directory = 'outputs/meta'

    with open('{}/dvc.md'.format(emitter_directory), 'w') as f:
        f.write(HUGO_HEADER)
        f.write('# Data Sources Summary\n')
        #f.write(tabulate(papers, headers="keys", tablefmt="pipe"))
        f.write(cardify(papers))
        f.write('\n')

        path = '{}/source*.dvc'.format(source_dir)

        source_data = []

        for filename in glob.iglob(path, recursive=True):
            with open(filename, 'r') as stream:
                dvc = yaml.load(stream)
                if 'cmd' not in dvc:
                    dvc['cmd'] = ''
                dvc = types.SimpleNamespace(**dvc)
                if hasattr(dvc, 'outs'):
                    for out in dvc.outs:
                        del out['cache']
                        out['path'] = out['path'].replace('_', '\\_')
                        source_data.append(out)

        f.write('# Data Sources Detail\n')
        f.write(tabulate(source_data, headers="keys", tablefmt="pipe"))
        f.write('\n')

        path = '{}/output*.dvc'.format(source_dir)

        out_data = []

        for filename in glob.iglob(path, recursive=True):
            with open(filename, 'r') as stream:
                dvc = yaml.load(stream)
                if 'cmd' not in dvc:
                    dvc['cmd'] = ''
                dvc = types.SimpleNamespace(**dvc)
                if hasattr(dvc, 'outs'):
                    for out in dvc.outs:
                        out['path'] = out['path'].replace('_', '\\_')
                    out_data.append({'cmd': dvc.cmd.replace('_', '\\_').split(';')[0], 'output': '<br/>'.join([(out['path']) for out in dvc.outs])})

        f.write('# Transformation Detail\n')
        f.write(tabulate(out_data, headers="keys", tablefmt="pipe"))
        f.write('\n')


if __name__ == '__main__':
    transform()

