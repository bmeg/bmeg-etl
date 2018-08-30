import pytest
import os
import json


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
        with open(vertex_file_path, 'r', encoding='utf-8') as f:
            for line in f:
                # should be json
                vertex_dict = json.loads(line)
                # should have all vertex keys
                Helpers.assert_vertex_keys_populated(vertex_dict)
                Helpers.assert_data_keys_populated(data_class, vertex_dict)

    @staticmethod
    def assert_edge_file_valid(from_data_class, to_data_class, edge_file_path):
        """ ensure file exists; populated with json objs with required keys """
        if not isinstance(from_data_class, list):
            from_data_class = [from_data_class]
        error_message = 'data_class {} {} {}' \
                        .format(from_data_class, to_data_class, edge_file_path)
        assert os.path.isfile(edge_file_path), error_message
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

    @staticmethod
    def assert_edge_joins_valid(graph_file_paths, exclude_labels=[]):
        """ load an in memory 'graph', ensure all edges link to a vertex and vice versa"""
        vertices = {}
        edges = {}
        for graph_file_path in graph_file_paths:
            with open(graph_file_path, 'r', encoding='utf-8') as f:
                store = vertices
                if 'Edge' in graph_file_path:
                    store = edges
                for line in f:
                    obj = json.loads(line)
                    store[obj['gid']] = obj
        # ensure that all edges have vertexes
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
            assert vertices[_from], 'edge {} from {} does not exist'.format(edge_gid, _from)
            assert vertices[_to], 'edge {} from {} does not exist'.format(edge_gid, _to)
        # ensure that all vertexes have edge
        froms = [edges[gid]['from'] for gid in edges.keys()]
        tos = [edges[gid]['to'] for gid in edges.keys()]
        for vertex_gid in vertices.keys():
            label = vertex_gid.split(':')[0]
            if label == 'Deadletter':
                continue
            assert vertex_gid in froms or vertex_gid in tos, 'could not find {} in edge'.format(vertex_gid)


@pytest.fixture
def helpers():
    return Helpers