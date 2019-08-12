import json
from transform.gdc.gdc_transformer import make_reader
from attrdict import AttrDict
from transform.pdc.pdc_harvester import load
from transform.tcia.tcia_harvester import load  as tcia_loader
import networkx as nx
from collections import defaultdict


def subgraph(node_pairs, graph):
    """Create a subgraph of the node and its pedigree."""
    nodes = []
    for node_pair in node_pairs:
        nodes.extend(node_pair)
    pedigree_nodes = []
    for node in set(nodes):
        pedigree_nodes.extend(
            [node] + [n for n in nx.ancestors(graph, node)] + [n for n in nx.descendants(graph, node)]
        )
    return graph.subgraph(pedigree_nodes)


class CDATransformer():

    def __init__(self, tcia, pdc, gdc):
        """Initializes self from passed graphs."""
        self.tcia = tcia
        self.pdc = pdc
        self.gdc = gdc

    def graph(self):
        """Adds `same_as` edge on `Case` vertex in composed cga networkx graph."""

        g = nx.compose(self.tcia,  nx.compose(self.pdc, self.gdc))

        # load cases
        gdc_cases = [ AttrDict(json.loads(c)) for c in make_reader('Case.Vertex')]
        pdc_cases = load('allCases')
        tcia_cases = [o for o in tcia_loader('patients')]

        # create maps
        gdc_submitter_id_map =  {
            c.data.submitter_id: {
                'gid': c.gid,
                'id': c.gid.replace('Case:','')
            } for c in gdc_cases
        }

        gdc_id_map =  {
            c.gid.replace('Case:',''): {
                'gid': c.gid,
                'submitter_id': c.data.submitter_id,
            } for c in gdc_cases
        }


        pdc_case_ids = set([c.gdc_case_id for c in pdc_cases])
        pdc_id_map = {
            c.gdc_case_id: c.case_submitter_id
            for c in pdc_cases
        }

        tcia_case_ids = set([c.PatientID for c in tcia_cases])

        # calculate intersections
        pdc_gdc_intersection = pdc_id_map.keys() & gdc_id_map.keys()
        tcia_gdc_intersection = tcia_case_ids & gdc_submitter_id_map.keys()
        pdc_tcia_gdc_intersection = tcia_case_ids & gdc_submitter_id_map.keys() & pdc_id_map.keys()

        print(
            'pdc_gdc_intersection:', len(pdc_gdc_intersection),
            'tcia_gdc_intersection:', len(tcia_gdc_intersection),
            'pdc_tcia_gdc_intersection:', len(pdc_tcia_gdc_intersection)
        )

        # add intersection to graph
        for k in pdc_gdc_intersection:
            pdc_id = pdc_id_map[k]
            gdc_id = gdc_id_map[k]
            g.add_edge(gdc_id['gid'], pdc_id, label='same_as')

        for k in tcia_gdc_intersection:
            tcia_id = k
            gdc_id = gdc_submitter_id_map[k]
            g.add_edge(gdc_id['gid'], tcia_id, label='same_as')

        self.cda = g
        return self.cda

    def shared_cases_summary(self):
        """Summarizes of shared_cases."""
        same_as_cases = set([(u,v) for u,v,d in self.cda.edges.data() if d['label'] == 'same_as'])

        def precedents(_id):
            path = 'cases.projects'.split('.')
            r = {}
            for label in path:
                for source, target, data in self.cda.in_edges(_id, data=True):
                    if data['label'] == label:
                        r[self.cda.node[_id]['label']] = _id
                        _id = source
            r[self.cda.node[_id]['label']] = _id
            return r

        project_shared_counts = defaultdict(int)
        for same_as in same_as_cases:
            precedents_0 = precedents(same_as[0])
            precedents_1 = precedents(same_as[1])
            key_0 = '{}.{}'.format(precedents_0['Source'], precedents_0['Project'])
            key_1 = '{}.{}'.format(precedents_1['Source'], precedents_1['Project'])
            key = '{}/{}'.format(*sorted([key_0, key_1]))
            project_shared_counts[key] += 1

        return project_shared_counts

    def same_as_graph(self):
        """Constructs a subgraph with only 'same_as' cases"""
        same_as_cases = set([(u,v) for u,v,d in self.cda.edges.data() if d['label'] == 'same_as'])
        return subgraph(same_as_cases, self.cda)
