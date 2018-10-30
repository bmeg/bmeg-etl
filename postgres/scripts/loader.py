<<<<<<< HEAD
<<<<<<< HEAD
=======
>>>>>>> sync
import logging
import yaml
import os
import types
import argparse

"""
Starts producers that read vertex and edge json files and write to FIFO queues
Starts consumers who read from FIFO queues and write to postgres tables using psql's COPY
"""


def start_consumer(consumer, psql_cmd):
    """ read from fifo, write to postgres, yields commands """
    consumer = types.SimpleNamespace(**consumer)
    yield '{} -c "{}"'.format(psql_cmd, consumer.cmd.strip().replace('"', '\\"'))


def transform_and_load(command, config):
    """ read from file, transform, write to fifo, yields commands """
    command = types.SimpleNamespace(**command)
    psql_cmd = config.psql_cmd
    for file in getattr(config, command.files):
        if not os.path.isfile(file):
            logging.warning('{} does not exist'.format(file))
            continue
        cat = 'cat'
        if 'gz' in file:
            cat = 'zcat'
        consumer = '{} -c "{}"'.format(psql_cmd, command.consumer.strip().replace('"', '\\"'))
        producer = command.producer.strip()
        yield '{} {} | {} | {}'.format(cat, file, producer, consumer)


def execute(ddls, psql_cmd):
    """ for each ddl write it to psql, yields commands """
    for ddl in ddls:
        yield 'cat << EOF | {}\n{}\nEOF\n'.format(psql_cmd, ddl)


"""
loads json vertex files into postgres table edge(gid, label, from, to, data)
"""


def main(drop, index, config):

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
    transform_and_load_commands = []
    index_commands = []
    if drop:
        setup_commands.extend(execute([config.drop], config.psql_cmd))
    else:
        logging.debug('skipping dropping tables')

    setup_commands.extend(execute([config.ddl], config.psql_cmd))

    for c in config.commands:
        transform_and_load_commands.extend(transform_and_load(c, config))

    if index:
        index_commands.extend(execute(config.indexes, config.psql_cmd))

    outputs = {
        'setup_commands.txt': setup_commands,
        'transform_and_load_commands.txt': transform_and_load_commands,
    }
    for path, commands in outputs.items():
        with open(path, 'w') as outfile:
            outfile.write("\n".join(commands))
            outfile.write("\n")
            logging.info('wrote {}'.format(path))

    c = 0
    for command in index_commands:
        path = 'index_{}.txt'.format(c)
        c += 1
        with open(path, 'w') as outfile:
            outfile.write(command)
            outfile.write("\n")
            logging.info('wrote {}'.format(path))
    path = 'index_commands.txt'
    with open(path, 'w') as outfile:
        c = 0
        for command in index_commands:
            outfile.write("bash index_{}.txt\n".format(c))
            c += 1
        logging.info('wrote {}'.format(path))

    logging.info("""
    Done.  To import:
    ```
    bash setup_commands.txt
    parallel --jobs 5 < transform_and_load_commands.txt
    parallel --jobs 5 < index_commands.txt
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
    args = parser.parse_args()
    logging.debug(vars(args))
    main(**vars(args))
<<<<<<< HEAD
=======
import logging
import yaml
import os
import types
import argparse

"""
Starts producers that read vertex and edge json files and write to FIFO queues
Starts consumers who read from FIFO queues and write to postgres tables using psql's COPY
"""


def start_consumer(consumer, psql_cmd):
    """ read from fifo, write to postgres, yields commands """
    consumer = types.SimpleNamespace(**consumer)
    yield '{} -c "{}"'.format(psql_cmd, consumer.cmd.strip().replace('"', '\\"'))


def transform_and_load(command, config):
    """ read from file, transform, write to fifo, yields commands """
    command = types.SimpleNamespace(**command)
    psql_cmd = config.psql_cmd
    for file in getattr(config, command.files):
        if not os.path.isfile(file):
            logging.warning('{} does not exist'.format(file))
            continue
        cat = 'cat'
        if 'gz' in file:
            cat = 'zcat'
        consumer = '{} -c "{}"'.format(psql_cmd, command.consumer.strip().replace('"', '\\"'))
        producer = command.producer.strip()
        yield '{} {} | {} | {}'.format(cat, file, producer, consumer)


def execute(ddls, psql_cmd):
    """ for each ddl write it to psql, yields commands """
    for ddl in ddls:
        yield 'cat << EOF | {}\n{}\nEOF\n'.format(psql_cmd, ddl)


"""
loads json vertex files into postgres table edge(gid, label, from, to, data)
"""


def main(drop, index, config):

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
    transform_and_load_commands = []
    index_commands = []
    if drop:
        setup_commands.extend(execute([config.drop], config.psql_cmd))
    else:
        logging.debug('skipping dropping tables')

    setup_commands.extend(execute([config.ddl], config.psql_cmd))

    for c in config.commands:
        transform_and_load_commands.extend(transform_and_load(c, config))

    if index:
        index_commands.extend(execute(config.indexes, config.psql_cmd))

    outputs = {
        'setup_commands.txt': setup_commands,
        'transform_and_load_commands.txt': transform_and_load_commands,
    }
    for path, commands in outputs.items():
        with open(path, 'w') as outfile:
            outfile.write("\n".join(commands))
            outfile.write("\n")
            logging.info('wrote {}'.format(path))

    c = 0
    for command in index_commands:
        path = 'index_{}.txt'.format(c)
        c += 1
        with open(path, 'w') as outfile:
            outfile.write(command)
            outfile.write("\n")
            logging.info('wrote {}'.format(path))
    path = 'index_commands.txt'
    with open(path, 'w') as outfile:
        c = 0
        for command in index_commands:
            outfile.write("bash index_{}.txt\n".format(c))
            c += 1
        logging.info('wrote {}'.format(path))

    logging.info("""
    Done.  To import, as a database user, run these commands :
    ```
    bash setup_commands.txt
    parallel --jobs 10 < transform_and_load_commands.txt
    parallel --jobs 10 < index_commands.txt
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
    args = parser.parse_args()
    logging.debug(vars(args))
    main(**vars(args))
>>>>>>> initial
=======
>>>>>>> sync
