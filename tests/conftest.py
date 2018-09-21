import pytest
import os
import json
import gripql
import contextlib


class Helpers:

    @staticmethod
    def assert_data_keys_populated(data_class, vertex_dict):
        """ ensure that all non Union(NoneType,...) fields are not empty. """
        # introspect mandatory keys
        for k in data_class.__dataclass_fields__.keys():
            field = data_class.__dataclass_fields__[k]
            # skip if union(None, ...)
            if 'typing.Union' in str(field.type) and 'NoneType' in str(field.type):
                continue
            assert vertex_dict['data'][k], 'empty key %s' % k

    @staticmethod
    def assert_vertex_keys_populated(vertex_dict):
        """ ensure that graph keys populated """
        # minimum graph keys
        assert list(vertex_dict.keys()) == ['_id', 'gid', 'label', 'data'], \
            'expected keys'
        for k in vertex_dict.keys():
            assert vertex_dict[k], 'empty key %s' % k

    @staticmethod
    def assert_vertex_file_valid(data_class, vertex_file_path):
        """ ensure file exists; populated with json objs with required keys """
        error_message = 'data_class {} {} ' \
                        .format(data_class, vertex_file_path)
        assert os.path.isfile(vertex_file_path), error_message
        c = 0
        with open(vertex_file_path, 'r', encoding='utf-8') as f:
            for line in f:
                # should be json
                vertex_dict = json.loads(line)
                # should have all vertex keys
                Helpers.assert_vertex_keys_populated(vertex_dict)
                Helpers.assert_data_keys_populated(data_class, vertex_dict)
                c += 1
        return c

    @staticmethod
    def assert_edge_file_valid(from_data_class, to_data_class, edge_file_path):
        """ ensure file exists; populated with json objs with required keys """
        if not isinstance(from_data_class, list):
            from_data_class = [from_data_class]
        error_message = 'data_class {} {} {}' \
                        .format(from_data_class, to_data_class, edge_file_path)
        assert os.path.isfile(edge_file_path), error_message
        c = 0
        with open(edge_file_path, 'r', encoding='utf-8') as f:
            for line in f:
                # should be json
                edge_dict = json.loads(line)
                # from, to should exist with correct gid
                found = False
                for clazz in from_data_class:
                    name = str(clazz.__name__).split('\\.')[-1]
                    if name in edge_dict['from']:
                        found = True
                assert found, 'edge.from should contain {} {}'.format(from_data_class, edge_dict['from'])

                name = str(to_data_class.__name__).split('\\.')[-1]
                assert name in edge_dict['to'], 'edge.to should contain {}'.format(name)
                c += 1
        return c

    @staticmethod
    def load_stores(graph_file_paths):
        """ load an in memory 'graph' returns (vertices, edges)"""
        vertices = {}
        edges = {}
        for graph_file_path in graph_file_paths:
            with contextlib.suppress(FileNotFoundError):
                with open(graph_file_path, 'r', encoding='utf-8') as f:
                    store = vertices
                    if 'Edge' in graph_file_path:
                        store = edges
                    for line in f:
                        obj = json.loads(line)
                        store[obj['gid']] = obj
        return vertices, edges

    @staticmethod
    def assert_edge_has_vertex(vertices, edges, exclude_labels=[]):
        """ensure that all edges have vertexes"""
        for edge_gid in edges.keys():
            edge = edges[edge_gid]
            _from = edge['from']  # from keyword workaround
            label = _from.split(':')[0]
            if label in exclude_labels:
                continue
            _to = edge['to']
            label = _to.split(':')[0]
            if label in exclude_labels:
                continue
            assert vertices.get(_from, None), 'edge {} from {} does not exist'.format(edge_gid, _from)
            assert vertices.get(_to, None), 'edge {} from {} does not exist'.format(edge_gid, _to)

    @staticmethod
    def assert_vertex_has_edge(vertices, edges, exclude_labels=[]):
        """ensure that all vertexes have edge"""
        froms = [edges[gid]['from'] for gid in edges.keys()]
        tos = [edges[gid]['to'] for gid in edges.keys()]
        for vertex_gid in vertices.keys():
            label = vertex_gid.split(':')[0]
            if label == 'Deadletter':
                continue
            assert vertex_gid in froms or vertex_gid in tos, 'could not find {} in edge'.format(vertex_gid)

    @staticmethod
    def assert_edge_joins_valid(graph_file_paths, exclude_labels=[]):
        """ load an in memory 'graph', ensure all edges link to a vertex and vice versa"""
        vertices, edges = Helpers.load_stores(graph_file_paths)
        Helpers.assert_edge_has_vertex(vertices, edges, exclude_labels)
        Helpers.assert_vertex_has_edge(vertices, edges, exclude_labels)


@pytest.fixture
def helpers():
    """ ratify Helper class """
    return Helpers


@pytest.fixture(scope="module")
def graph():
    """ return a connection to the bmeg graph (control w/ BMEG_URL, BMEG_GRAPH env var) """
    bmeg_url = os.getenv('BMEG_URL', "http://arachne.compbio.ohsu.edu")
    bmeg_graph = os.getenv('BMEG_GRAPH', "bmeg_test")
    return gripql.Connection(bmeg_url).graph(bmeg_graph)


@pytest.fixture(scope="module")
def V(graph):
    """ return the Vertex in the bmeg graph """
    return graph.query().V()


@pytest.fixture
def emitter_directory(request):
    """ get the full path of the test output """
    return os.path.join(request.fspath.dirname, 'test')


@pytest.fixture
def emitter_prefix(request):
    """ get the full path of the test output """
    return 'test'
