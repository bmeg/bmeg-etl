import logging
import yaml
import os
import types
import argparse
from jinja2 import Template

"""
creates a set of commands for loading and indexing graph
"""


def edge_labels(file):
    """ given an Edge file path, get it's from and to """
    LABELS = {
        'AliquotFor': {'from': 'Aliquot', 'to': 'Biosample'},
        'AlleleCall': {'from': 'Allele', 'to': 'Callset'},
        'AlleleIn': {'from': 'Allele', 'to': 'Gene'},
        'BiosampleFor': {'from': 'Biosample', 'to': 'Individual'},
        'COCAClusterFor': {'from': 'COCACluster', 'to': 'Individual'},
        'CallsetFor': {'from': 'Callset', 'to': 'Aliquot'},
        'DrugResponseIn': {'from': 'DrugResponse', 'to': 'Aliquot'},
        'ExonFor': {'from': 'Exon', 'to': 'Transcript'},
        'ExpressionOf': {'from': 'Expression', 'to': 'Aliquot'},
        'GeneOntologyAnnotation': {'from': 'GO', 'to': 'PFAM'},
        'GeneOntologyIsA': {'from': 'GO', 'to': 'GO'},
        'HasAlleleFeature': {'from': 'G2PAssociation', 'to': 'Allele'},
        'HasEnvironment': {'from': 'G2PAssociation', 'to': 'Compound'},
        'HasGeneFeature': {'from': 'G2PAssociation', 'to': 'Gene'},
        'HasMinimalAlleleFeature': {'from': 'G2PAssociation', 'to': 'MinimalAllele'},
        'HasPhenotype': {'from': 'G2PAssociation', 'to': 'Phenotype'},
        'HasSupportingReference': {'from': 'G2PAssociation', 'to': 'Publication'},
        'InProject': {'from': 'Individual', 'to': 'Project'},
        'MinimalAlleleIn': {'from': 'MinimalAllele', 'to': 'Gene'},
        'PFAMAlignment': {'from': 'Protein', 'to': 'PFAM'},
        'PFAMClanMember': {'from': 'PFAMCLAN', 'to': 'PFAM'},
        'PhenotypeOf': {'from': 'Aliquot', 'to': 'Phenotype'},
        'ProteinFor': {'from': 'Protein', 'to': 'Transcript'},
        'ResponseTo': {'from': 'DrugResponse', 'to': 'Compound'},
        'StructureFor': {'from': 'PDB', 'to': 'Protein'},
        'TranscriptFor': {'from': 'Transcript', 'to': 'Gene'},
        'TreatedWith': {'from': 'Individual', 'to': 'Compound'}
    }

    if 'Edge' not in file:
        return None, None

    file_parts = '.'.join(file.split('/')).split('.')
    label_index = file_parts.index('Edge') - 1
    label = file_parts[label_index]
    return ':{}'.format(LABELS[label]['from']), ':{}'.format(LABELS[label]['to'])


def execute(ddls, cypher_shell):
    """ for each ddl write it to cypher_shell, yields commands """
    for ddl in ddls:
        yield 'cat << EOF | {}\n{}\nEOF\n'.format(cypher_shell, ddl)


def transform_and_load(command, template, files, limit):
    """ read from file, transform, write to fifo, yields commands """
    t = Template(template)
    if limit:
        limit = 'LIMIT {}'.format(limit)
    else:
        limit = ''
    for file in files:
        if not os.path.isfile(file):
            logging.warning('{} does not exist'.format(file))
            continue
        from_label, to_label = edge_labels(file)
        url = 'file:///{}'.format(file)
        for cmd in execute([t.render(url=url, limit=limit, from_label=from_label, to_label=to_label)], command):
            yield cmd


"""
loads json vertex files into postgres table edge(gid, label, from, to, data)
"""


def main(drop, index, config, limit):

    with open(config, 'r') as stream:
        config = yaml.load(stream)

    config = types.SimpleNamespace(**config)
    config.edge_files = config.edge_files.strip().split()
    config.vertex_files = config.vertex_files.strip().split()
    config.matrix_files = config.matrix_files.strip().split()

    # # ensure files exist
    # all_files = config.vertex_files + config.edge_files + config.matrix_files
    # logging.debug(all_files)
    # for fname in all_files:
    #     assert os.path.isfile(fname), '{} does not exist'.format(fname)

    setup_commands = []
    vertex_commands = []
    index_commands = []
    edge_commands = []

    if drop:
        setup_commands.extend(execute([config.drop], config.cypher_shell))
    else:
        logging.debug('skipping dropping tables')

    setup_commands.extend(execute([config.ddl], config.cypher_shell))

    vertex_commands.extend(transform_and_load(command=config.cypher_shell, template=config.vertex_template, files=config.vertex_files, limit=limit))
    edge_commands.extend(transform_and_load(command=config.cypher_shell, template=config.edge_template, files=config.edge_files, limit=limit))

    if index:
        index_commands.extend(execute(config.indexes, config.cypher_shell))

    path = 'setup_commands.txt'
    with open(path, 'w') as outfile:
        outfile.write("\n".join(setup_commands))
        outfile.write("\n")
        logging.info('wrote {}'.format(path))

    c = 0
    for command in vertex_commands:
        path = 'vertex_{}.txt'.format(c)
        c += 1
        with open(path, 'w') as outfile:
            outfile.write(command)
            outfile.write("\n")

    path = 'vertex_commands.txt'
    with open(path, 'w') as outfile:
        c = 0
        for command in vertex_commands:
            outfile.write("bash vertex_{}.txt\n".format(c))
            c += 1
        logging.info('wrote {}'.format(path))

    c = 0
    for command in index_commands:
        path = 'index_{}.txt'.format(c)
        c += 1
        with open(path, 'w') as outfile:
            outfile.write(command)
            outfile.write("\n")

    path = 'index_commands.txt'
    with open(path, 'w') as outfile:
        c = 0
        for command in index_commands:
            outfile.write("bash index_{}.txt\n".format(c))
            c += 1
        logging.info('wrote {}'.format(path))

    c = 0
    for command in edge_commands:
        path = 'edge_{}.txt'.format(c)
        c += 1
        with open(path, 'w') as outfile:
            outfile.write(command)
            outfile.write("\n")

    path = 'edge_commands.txt'
    with open(path, 'w') as outfile:
        c = 0
        for command in edge_commands:
            outfile.write("bash edge_{}.txt\n".format(c))
            c += 1
        logging.info('wrote {}'.format(path))

    logging.info("""
    Done.  To import:
    ```
    bash setup_commands.txt
    parallel --jobs 5 < vertex_commands.txt
    parallel --jobs 5 < index_commands.txt
    parallel --jobs 5 < edge_commands.txt
    ```
    """)


if __name__ == '__main__':  # pragma: no cover
    # log setup
    logging.getLogger().setLevel(logging.DEBUG)

    parser = argparse.ArgumentParser(description="Loads vertexes and edges into postgres")
    parser.add_argument('--drop', dest='drop', action='store_true', default=False, help="drop the tables first [False]")
    parser.add_argument('--skip_index', dest='index', action='store_false', default=True, help="index after loading [True]")
    config_path = "{}/config.yml".format(os.path.dirname(os.path.realpath(__file__)))
    parser.add_argument('--config', dest='config', default=config_path, help="config path {}".format(config_path))
    parser.add_argument('--limit', dest='limit', type=int, default=None, help="limit the number of rows in each vertex/edge")
    args = parser.parse_args()
    logging.debug(vars(args))
    main(**vars(args))
