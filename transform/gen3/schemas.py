from bmeg.ioutils import reader
import bmeg.vertex
from dataclasses import fields

import sys
import os
import json
from collections import defaultdict
import subprocess

import networkx as nx
import matplotlib.pyplot as plt

import inflection

from yaml import dump, load, FullLoader, Dumper

import copy

if sys.version > '3':
    long = int


# https://github.com/yaml/pyyaml/issues/103
class NoAliasDumper(Dumper):
    def ignore_aliases(self, data):
        return True


def gid_exceptions(label):
    """Known exceptions to Label:Id convention."""
    if label.startswith('ENSG'):
        return 'Gene'
    if label.startswith('ENSP'):
        return 'Protein'
    if label.startswith('ENST'):
        return 'Transcript'
    if label.startswith('ENSE'):
        return 'Exon'
    if label.startswith('PDB'):
        return 'ProteinStructure'
    if label.startswith('GO'):
        return 'GeneOntologyTerm'
    return label


def edge_record():
    """Stores edge info"""
    def record():
        return {'src': None, 'dst': None, 'label': None, 'files': []}
    return defaultdict(record)


def deduce_edges(find_cmd, output_file):
    """Samples all edges in source.
    stores:
        {
          "<label>": {
            "<dst>": {
              "src": "<vertex>",
              "dst": "<vertex>",
              "label": "<vertex>",
              "files": [
                "<edge-path>",
                "<edge-path>",
                "<edge-path>",
    """
    p = subprocess.Popen(find_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    p.wait()

    edges = defaultdict(edge_record)
    for line in p.stdout.readlines():
        line = line.rstrip().decode('UTF-8')
        if not os.path.isfile(line):
            print('# WARNING: ', line, 'is not a file')
            continue
        with reader(line) as ins:
            for edge in ins:
                edge = json.loads(edge)
                src = gid_exceptions(edge['from'].split(':')[0])
                dst = gid_exceptions(edge['to'].split(':')[0])
                label = edge['label']
                edges[label][dst]['src'] = src
                edges[label][dst]['dst'] = dst
                edges[label][dst]['label'] = label
                edges[label][dst]['files'].append(line)
                break

    with open(output_file, 'w') as output:
        output.write(json.dumps(edges))
    return output_file


def load_edges(file_path):
    """Loads the edges file"""
    return json.load(open(file_path, 'r'))


def create_graph(edges):
    """Creates a graph from all edges"""
    G = nx.MultiDiGraph()
    nodes = set([])
    for lable in edges:
        for k, edge in edges[lable].items():
            nodes.add(edge['src'])
            nodes.add(edge['dst'])
    G.add_nodes_from(nodes)
    for lable in edges:
        for k, edge in edges[lable].items():
            G.add_edge(edge['src'], edge['dst'], **edge)
    return G


def draw_graph(G, file_name='bmeg.png'):
    """Creates a shell plot."""
    pos = nx.shell_layout(G)
    nx.draw(G, pos, node_color='#A0CBE2', font_size=8,
            width=1, edge_cmap=plt.cm.Blues, with_labels=True)
    plt.savefig(file_name)
    return file_name


def create_schema(G):
    """Compares G to bmeg.vertex"""
    bmeg_vertices = set([x for x in dir(bmeg.vertex) if x[0].isupper()])
    deduced_vertices = set([n for n in G.nodes])
    print('# INFO bmeg vertex without vertex deduced from edge files.', bmeg_vertices - deduced_vertices)
    print('# INFO deduced vertices without bmeg vertex', deduced_vertices - bmeg_vertices)
    print('# INFO vertices in common', deduced_vertices & bmeg_vertices)

    def schema_type(python_type):
        type_map = {
            bool: 'boolean',
            long: 'integer',
            float: 'number',
            dict: 'object',
            list: 'array',
            str: 'string',
            type(None): 'null'
        }
        if python_type in type_map:
            return {'type': type_map.get(python_type)}
        if 'typing.Union' in str(python_type):
            return {'type': [type_map[a] for a in python_type.__args__]}
        if 'enum' in str(python_type):
            return {'enum': [k.name for k in python_type]}
        print('# WARNING: ', 'no type for', python_type)
        return '{}'

    schema = {}
    for v in deduced_vertices & bmeg_vertices:
        field_types = {field.name: schema_type(field.type) for field in fields(eval('bmeg.vertex.{}'.format(v)))}
        schema[v] = {}
        schema[v]['links'] = []
        schema[v]['properties'] = field_types
        schema[v]['required'] = [k for k, o in schema[v]['properties'].items() if type(o['type']) is not list]

    to_many = {'type': {'$ref': '_definitions.yaml#/to_many'}}
    for v in deduced_vertices & bmeg_vertices:
        source_edges = [e for e in G.edges if e[0] == v]
        # [('Case', 'Project', 0), ('Case', 'Compound', 0)]
        for src, dst, junk in source_edges:
            edge_data = G.get_edge_data(src, dst)
            for i, e in edge_data.items():
                # e {'src': 'Case', 'dst': 'Compound', 'label': 'TreatedWith', 'files': ['outputs/compound/normalized.TreatedWith.Edge.json.gz', 'outputs/gdc/TreatedWith.Edge.json.gz']}
                schema[v]['links'].append(
                    {
                        "backref": inflection.underscore(inflection.pluralize(e['src'])),
                        "label": inflection.underscore(e['label']),
                        "multiplicity": "many_to_many",
                        "name": inflection.underscore(inflection.pluralize(e['dst'])),
                        "required": False,
                        "target_type": inflection.underscore(e['dst'])
                    }
                )
                schema[v]['properties'][inflection.underscore(inflection.pluralize(e['dst']))] = to_many
                if e['dst'] in schema:
                    schema[e['dst']]['properties'][inflection.underscore(inflection.pluralize(e['src']))] = to_many
    return schema


def create_yaml(schema, output_dir='outputs/gen3'):
    """Creates yaml files."""
    this_dir = os.path.dirname(os.path.realpath(__file__))
    default = {}
    with open(os.path.join(this_dir, 'default.yaml')) as ins:
        default = load(ins, Loader=FullLoader)
    paths = []
    for label, vertex in schema.items():
        schema = copy.deepcopy(default)
        schema['properties'].update(vertex['properties'])
        schema['links'] = schema['links'] + vertex['links']
        schema['required'] = list(set(schema['required'] + vertex['required']))
        schema_id = inflection.underscore(label)
        schema['title'] = label
        schema['description'] = 'Autogenerated definitions for {}'.format(schema['title'])
        schema['id'] = schema_id
        file_name = '{}.yaml'.format(schema_id)
        p = os.path.join(output_dir, file_name)
        with open(p, 'w') as outs:
            outs.write(dump(schema, default_flow_style=False, Dumper=NoAliasDumper))
        paths.append(p)
    return paths


def transform(find_cmd='cat scripts/bmeg_file_manifest.txt | grep Edge', file_path='source/gen3/edges.json', output_dir='outputs/gen3'):
    """Creates schema files."""
    G = create_graph(load_edges(deduce_edges(find_cmd, file_path)))
    schema = create_schema(G)
    return create_yaml(schema, output_dir)


if __name__ == "__main__":
    transform()
