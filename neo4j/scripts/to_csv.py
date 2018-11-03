import ujson
import logging
import yaml
import types
import os
import csv
from bmeg.ioutils import reader
import argparse


CONTROL = dict(zip(range(32), ' ' * 32))


def strip_control(s):
    """ clean string """
    return s.translate(CONTROL)


def to_vertex(path):
    """ vertex with only scalar data"""
    with reader(path) as ins:
        for line in ins:
            line = ujson.loads(line)
            del line['_id']
            for k, v in line['data'].items():
                if isinstance(v, dict) or isinstance(v, list):
                    continue
                if isinstance(line['data'][k], str):
                    line[k] = strip_control(line['data'][k])
                else:
                    line[k] = line['data'][k]
            # remove data map and label
            del line['data']
            del line['label']
            if 'source_document' in line:
                del line['source_document']
            yield line


def to_edge(path):
    """ edge with only scalar data"""
    with reader(path) as ins:
        for line in ins:
            line = ujson.loads(line)
            # not used
            del line['_id']
            del line['gid']
            # move data from map to root
            for k, v in line['data'].items():
                if isinstance(v, dict) or isinstance(v, list):
                    continue
                if isinstance(line['data'][k], str):
                    line[k] = strip_control(line['data'][k])
                else:
                    line[k] = line['data'][k]
            # remove data map
            del line['data']
            # rename from, to, label
            line[':START_ID'] = line['from']
            line[':END_ID'] = line['to']
            line[':TYPE'] = line['label']
            del line['label']
            del line['from']
            del line['to']
            yield line


def keys(path, sample_size=1000):
    """ return [names of keys] and [neo types of keys]"""
    if 'Expression' in path:
        sample_size = 1  # no need to read huge, uniform records
    keys = []
    values = []
    py_2_neo = {
        'str': 'string',
        'list': 'string[]'
    }  # xlate py types to neo4j
    c = 0
    xformer = to_vertex
    if 'Edge' in path:
        xformer = to_edge
    for line in xformer(path):
        if len(line.keys()) > len(keys):
            keys = line.keys()
        if len(line.values()) > len(values):
            values = line.values()
        if c == sample_size:
            break
        c += 1

    decorated_keys = []
    for key, value in zip(keys, values):
        py_type = value.__class__.__name__
        value_type = py_2_neo.get(py_type, py_type)
        if key == 'gid':
            value_type = 'ID'
        if key in [':START_ID', ':END_ID', ':TYPE']:
            decorated_keys.append(key)
        else:
            decorated_keys.append('{}:{}'.format(key, value_type))
    return keys, decorated_keys


def values(path):
    """ return a dict for each line """
    xformer = to_vertex
    if 'Edge' in path:
        xformer = to_edge
    for line in xformer(path):
        yield line


def get_label(path):
    """ given path, return lable """
    label_type = 'Vertex'
    if 'Edge' in path:
        label_type = 'Edge'
    file_parts = '.'.join(path.split('/')).split('.')
    label_index = file_parts.index(label_type) - 1
    return file_parts[label_index]


def to_csv_header(path, output=None):
    """ file to csv '{path}.headercsv' """
    fieldnames, decorated_fieldnames = keys(path)
    return dict(zip(fieldnames, decorated_fieldnames))


def get_output_path(path):
    return 'neo4j/scripts/{}.csv'.format(path.replace('/', '.'))


def to_csv(path, limit=None, header=False, output=None, header_dict=None):
    """ file to csv '{path}.csv' """
    if header_dict:
        fieldnames, decorated_fieldnames = header_dict.keys(), header_dict.values()
    else:
        fieldnames, decorated_fieldnames = keys(path)
    output_path = output
    if not output:
        output_path = 'neo4j/scripts/{}.header.csv'.format(path)
    with open(output_path, "w", newline='') as myfile:
        writer = csv.DictWriter(myfile, fieldnames=fieldnames, extrasaction='ignore')
        if header:
            writer.writeheader()
        c = 0
        for line in values(path):
            writer.writerow(line)
            c += 1
            if limit and c == limit:
                break
        logging.info(c)
    return output_path


def to_csv_job(path, limit=None):
    """ cmd line to transform json to csv """
    output_path = get_output_path(path)
    if limit:
        limit = '--limit {}'.format(limit)
    else:
        limit = ''
    return 'python neo4j/scripts/to_csv.py --input {} --output {} {}'.format(path, output_path, limit)


def main(config, limit, input, output):
    """ read config, render csv file and neo4j-import clause"""

    with open(config, 'r') as stream:
        config = yaml.load(stream)

    config = types.SimpleNamespace(**config)
    config.edge_files = list(set(config.edge_files.strip().split()))
    config.vertex_files = list(set(config.vertex_files.strip().split()))
    config.matrix_files = list(set(config.matrix_files.strip().split()))

    # input specified, just run and exit
    if input:
        header = {}
        label = get_label(input)
        if 'Edge' in input:
            files = config.edge_files
        else:
            files = config.vertex_files + config.matrix_files
        for path in files:
            if '.{}.'.format(label) not in path:
                continue
            header = {**header, **to_csv_header(path)}
        output_path = to_csv(path=input, output=output, limit=limit, header_dict=header)
        logging.info(output_path)
        return

    vertex_csvs = {}
    edge_csvs = {}
    headers = {}
    to_csv_commands = []
    # read all files to determine header by label
    for path in config.vertex_files + config.matrix_files + config.edge_files:
        if not os.path.isfile(path):
            logging.warning('{} does not exist'.format(path))
            continue
        label = get_label(path)
        if label not in headers:
            headers[label] = {}
        headers[label] = {**headers[label], **to_csv_header(path)}
    # write csv header files
    for label in headers.keys():
        output_path = 'neo4j/scripts/{}.header.csv'.format(label)
        with open(output_path, "w", newline='') as myfile:
            writer = csv.DictWriter(myfile, fieldnames=headers[label].keys())
            writer.writerow(headers[label])

    for path in config.vertex_files + config.matrix_files:
        if not os.path.isfile(path):
            continue
        label = get_label(path)
        if label not in vertex_csvs:
            vertex_csvs[label] = []
            vertex_csvs[label].append('neo4j/scripts/{}.header.csv'.format(label))
        to_csv_commands.append(to_csv_job(path, limit=limit))
        vertex_csvs[label].append(get_output_path(path))

    for path in config.edge_files:
        if not os.path.isfile(path):
            logging.warning('{} does not exist'.format(path))
            continue
        label = get_label(path)
        if label not in edge_csvs:
            edge_csvs[label] = []
            edge_csvs[label].append('neo4j/scripts/{}.header.csv'.format(label))
        to_csv_commands.append(to_csv_job(path, limit=limit))
        edge_csvs[label].append(get_output_path(path))

    path = 'neo4j/scripts/to_csv_commands.txt'
    with open(path, 'w') as outfile:
        for command in to_csv_commands:
            outfile.write("{}\n".format(command))
        logging.info('wrote {}'.format(path))

    nodes = []
    for key in vertex_csvs.keys():
        nodes.append('--nodes:{} {}'.format(key, ','.join(vertex_csvs[key])))

    edges = []
    for key in edge_csvs.keys():
        edges.append('--relationships:{} {}'.format(key, ','.join(edge_csvs[key])))

    cmds = '\n'.join([
        'parallel --jobs 10 < neo4j/scripts/to_csv_commands.txt',
        'sudo --user=neo4j neo4j-admin import --database test.db --ignore-missing-nodes=true --ignore-duplicate-nodes=true --ignore-extra-columns=true --high-io=true \\'
    ])
    logging.info('\n{}\n  {}\n'.format(cmds, ' \\\n  '.join(nodes + edges)))


if __name__ == '__main__':  # pragma: no cover
    logging.getLogger().setLevel(logging.DEBUG)
    parser = argparse.ArgumentParser(description="Loads vertexes and edges into postgres")
    config_path = "{}/config.yml".format(os.path.dirname(os.path.realpath(__file__)))
    parser.add_argument('--config', dest='config', default=config_path, help="config path {}".format(config_path))
    parser.add_argument('--limit', dest='limit', type=int, default=None, help="limit the number of rows in each vertex/edge")
    parser.add_argument('--input', dest='input', default=None, help="single input file")
    parser.add_argument('--output', dest='output', default=None, help="single output file")
    args = parser.parse_args()
    logging.debug(vars(args))
    main(**vars(args))
