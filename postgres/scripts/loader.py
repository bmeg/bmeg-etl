import logging
import yaml
import os
import types
import argparse

"""
Starts producers that read vertex and edge json files and write to FIFO queues
Starts consumers who read from FIFO queues and write to postgres tables using psql's COPY
"""


def ensure_fifo(path):
    """ create fifo, returns command """
    if not os.path.isfile(path):
        cmd = 'mkfifo {}'.format(path)
        return cmd


def start_consumer(consumer, psql_cmd):
    """ read from fifo, write to postgres, yields commands """
    consumer = types.SimpleNamespace(**consumer)
    yield ensure_fifo(consumer.fifo)
    yield '{} -c "{}" < {} &'.format(psql_cmd, consumer.cmd.strip().replace('"', '\\"'), consumer.fifo)


def start_producer(producer, config):
    """ read from file, transform, write to fifo, yields commands """
    producer = types.SimpleNamespace(**producer)
    for file in getattr(config, producer.files):
        if not os.path.isfile(file):
            logging.warning('{} does not exist'.format(file))
            continue
        cat = 'cat'
        if 'gz' in file:
            cat = 'zcat'
        yield '{} {} | head -3 | {} > {}'.format(cat, file, producer.cmd.strip(), producer.fifo)


def execute(ddls, psql_cmd):
    """ for each ddl write it to psql, yields commands """
    for ddl in ddls:
        yield 'cat << EOF | {}\n{}\nEOF\n'.format(psql_cmd, ddl)


"""
loads json vertex files into postgres table edge(gid, label, from, to, data)
"""


def main(dry, drop, index, config):

    with open(config, 'r') as stream:
        config = yaml.load(stream)

    config = types.SimpleNamespace(**config)
    config.edge_files = config.edge_files.strip().split()
    config.vertex_files = config.vertex_files.strip().split()
    config.matrix_files = config.matrix_files.strip().split()
    all_files = config.vertex_files + config.edge_files + config.matrix_files

    # # ensure files exist
    logging.debug(all_files)
    # for fname in all_files:
    #     assert os.path.isfile(fname), '{} does not exist'.format(fname)

    setup_commands = []
    consumer_commands = []
    producer_commands = []
    index_commands = []
    if drop:
        setup_commands.extend(execute([config.drop], config.psql_cmd))
    else:
        logging.debug('skipping dropping tables')

    setup_commands.extend(execute([config.ddl], config.psql_cmd))

    for c in config.consumers:
        consumer_commands.extend(start_consumer(c, config.psql_cmd))

    for p in config.producers:
        producer_commands.extend(start_producer(p, config))

    if index:
        index_commands.extend(execute([config.indexes], config.psql_cmd))

    outputs = {
        'setup_commands.txt': setup_commands,
        'consumer_commands.txt': consumer_commands,
        'producer_commands.txt': producer_commands,
        'index_commands.txt': index_commands
    }
    for path, commands in outputs.items():
        with open(path, 'w') as outfile:
            outfile.write("\n".join(commands))
            outfile.write("\n")
            logging.info('wrote {}'.format(path))

    logging.info('Done!')


if __name__ == '__main__':  # pragma: no cover
    # log setup
    logging.getLogger().setLevel(logging.DEBUG)

    parser = argparse.ArgumentParser(description="Loads vertexes and edges into postgres")
    parser.add_argument('--dry', dest='dry', action='store_true', default=False, help="echo commands, don't execute [False]")
    parser.add_argument('--drop', dest='drop', action='store_true', default=False, help="drop the tables first [False]")
    parser.add_argument('--skip_index', dest='index', action='store_false', default=True, help="index after loading [True]")
    config_path = "{}/config.yml".format(os.path.dirname(os.path.realpath(__file__)))
    parser.add_argument('--config', dest='config', default=config_path, help="config path {}".format(config_path))
    args = parser.parse_args()
    logging.debug(vars(args))
    main(**vars(args))
