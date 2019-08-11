from collections import defaultdict
from functools import reduce
import networkx as nx
import json
from bmeg.ioutils import reader
from attrdict import AttrDict

from bmeg.ioutils import reader
import json

def make_reader(name):
    """Instantiates a reader."""
    return reader('outputs/gdc/{}.json.gz'.format(name))

def vertex_readers():
    """Reader for all vertexes."""
#         Program.Vertex
    vertexes = """
    Aliquot.Vertex
    Case.Vertex
    Project.Vertex
    Sample.Vertex
    """.strip().split()    
    return [make_reader(v) for v in vertexes]
    
def edge_readers():
    """Reader for all edges."""
#     edges = """
#     Aliquot_Projects_Project.Edge
#     Aliquot_Sample_Sample.Edge
#     Case_Projects_Project.Edge
#     Case_Samples_Sample.Edge
#     Program_Projects_Project.Edge
#     Project_Aliquots_Aliquot.Edge
#     Project_Cases_Case.Edge
#     Project_Programs_Program.Edge
#     Project_Samples_Sample.Edge
#     Sample_Aliquots_Aliquot.Edge
#     Sample_Case_Case.Edge
#     Sample_Projects_Project.Edge    
#     """.strip().split() 
    edges = """
    Case_Samples_Sample.Edge
    Project_Cases_Case.Edge
    Sample_Aliquots_Aliquot.Edge
    """.strip().split()    
    return [make_reader(e) for e in edges]
    

class GDCTransformer():

    
    def graph(self):
        """Creates graph"""
        G = nx.MultiDiGraph()
        
        def add_nodes(g, vertexes):
            """Adds node to graph w/ label."""
            for v in vertexes:
                v = AttrDict(json.loads(v))
                g.add_node(v.gid, label=v.label)

        def add_edges(g, edges):
            """Adds edge to graph w/ label."""
            for e in edges:
                e = AttrDict(json.loads(e))
                g.add_edge(e['from'], e['to'], label=e.label)
                
                
        for vr in vertex_readers():
            add_nodes(G, vr)
            vr.close()

        for er in edge_readers():
            add_edges(G, er)
            er.close()
        
        case_reader = make_reader('Case.Vertex')
        for c in case_reader:
            c = AttrDict(json.loads(c))
            if c.data.gdc_attributes.diagnoses:
                for d in c.data.gdc_attributes.diagnoses:
                    G.add_node(d.diagnosis_id, label='Diagnosis')
                    G.add_edge(c.gid, d.diagnosis_id, label='diagnoses')
            d = c.data.gdc_attributes.demographic
            if d:
                G.add_node(d.demographic_id, label='Demographic')
                G.add_edge(c.gid, d.demographic_id, label='demographics')                
   
        G.add_node('gdc', label= 'Source')
        G.add_edges_from([('gdc',k) for (k, v) in G.nodes.data() if v['label'] == 'Project'], label='projects')
            
        return G
