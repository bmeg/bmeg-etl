from bmeg.ioutils import reader

import sys
import os
import json
from collections import defaultdict
import subprocess

import networkx as nx
import matplotlib.pyplot as plt

import matplotlib.pyplot
import pandas
import upsetplot


if sys.version > '3':
    long = int


CATEGORIES = {
    "acknowledgement": "administrative",
    "aligned_reads_index": "index_file",
    "aliquot": "biospecimen",
    "allele": "reference",
    "bcc_biomarker": "bcc_extention",
    "bcc_chemotherapy": "bcc_extention",
    "bcc_demographic": "bcc_extention",
    "bcc_diagnosis": "bcc_extention",
    "bcc_lesion": "bcc_extention",
    "bcc_participant": "bcc_extention",
    "bcc_radiotherapy": "bcc_extention",
    "bcc_sample": "bcc_extention",
    "bcc_surgery": "bcc_extention",
    "bcc_weight": "bcc_extention",
    "case": "administrative",
    "clinical_test": "clinical",
    "compound": "administrative",
    "core_metadata_collection": "administrative",
    "demographic": "clinical",
    "diagnosis": "clinical",
    "drug_response": "genomic",
    "experiment": "administrative",
    "experimental_metadata": "metadata_file",
    "exposure": "clinical",
    "family_history": "clinical",
    "gene": "reference",
    "gene_ontology_term": "reference",
    "genetrails_variant": "bcc_extention",
    "hop_survey": "hop extention",
    "keyword": "administrative",
    "observation": "clinical",
    "phenotype": "reference",
    "program": "administrative",
    "project": "administrative",
    "publication": "administrative",
    "read_group": "biospecimen",
    "read_group_qc": "notation",
    "sample": "biospecimen",
    "slide": "biospecimen",
    "slide_count": "notation",
    "slide_image": "data_file",
    "submitted_aligned_reads": "data_file",
    "submitted_copy_number": "data_file",
    "submitted_file": "data_file",
    "submitted_methylation": "data_file",
    "submitted_somatic_mutation": "data_file",
    "submitted_unaligned_reads": "data_file",
    "treatment": "clinical"
}


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
        if 'Edge' not in line:
            continue
        if not os.path.isfile(line):
            print('# WARNING: ', line, 'is not a file')
            continue
        with reader(line) as ins:
            c = 0
            for edge in ins:
                edge = json.loads(edge)
                c += 1
                src = gid_exceptions(edge['from'].split(':')[0])
                dst = gid_exceptions(edge['to'].split(':')[0])
                label = line
                edges[label][dst]['src'] = src
                edges[label][dst]['dst'] = dst
                edges[label][dst]['label'] = label
                edges[label][dst]['count'] = c
                # edges[label][dst]['files'].add(line)
            print(edges[label][dst])
    with open(output_file, 'w') as output:
        output.write(json.dumps(edges))
    return output_file


def count_vertexes(find_cmd, g):
    """Adds node count to graph"""
    p = subprocess.Popen(find_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    p.wait()
    vertex_counts = defaultdict(int)
    for line in p.stdout.readlines():
        line = line.rstrip().decode('UTF-8')
        if 'Vertex' not in line:
            continue
        if not os.path.isfile(line):
            print('# WARNING: ', line, 'is not a file')
            continue
        with reader(line) as ins:
            c = 0
            for vertex in ins:
                vertex = json.loads(vertex)
                c += 1
                vertex_counts[vertex['label']] = c
    vertex_counts = {k: f"{k}-{v}" for k, v in vertex_counts.items()}
    # ensure uncounted node has label
    for k in g.nodes():
        if k not in vertex_counts:
            vertex_counts[k] = f"{k}-?"
    g.vertex_counts = vertex_counts


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
    edge_labels = {}
    for lable in edges:
        for k, edge in edges[lable].items():
            print(lable, edge['src'], edge['dst'])
            G.add_edge(edge['src'], edge['dst'], **edge)
            edge_labels[(edge['src'], edge['dst'])] = edge['count']
    G.edge_labels = edge_labels
    return G


def draw_graph(G, file_name='bmeg.png'):
    """Creates a shell plot."""
    pos = nx.circular_layout(G)
    nx.draw(G, pos, node_color='#A0CBE2', font_size=8,
            width=1, edge_cmap=plt.cm.Blues, with_labels=True, labels=G.vertex_counts)
    nx.draw_networkx_edge_labels(G, pos, edge_labels=G.edge_labels)
    plt.savefig(file_name)
    return file_name


def draw_intersections():
    """Upset plots for genes and pubs."""
    p = subprocess.Popen("zcat outputs/pharmgkb/Publication_G2PAssociations_G2PAssociation.Edge.json.gz  | jq .from | sort -u", shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    pharmgkb_pubs = [line.strip() for line in p.stdout.readlines()]
    p = subprocess.Popen("zcat outputs/g2p/Publication_G2PAssociations_G2PAssociation.Edge.json.gz | jq .from | sort -u", shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    g2p_pubs = [line.strip() for line in p.stdout.readlines()]
    sample_df = pandas.DataFrame(upsetplot.from_contents({'pharmgkb': pharmgkb_pubs, 'g2p': g2p_pubs}))
    upsetplot.plot(sample_df, sort_by="cardinality", sum_over=False, show_counts='%d')
    current_figure = matplotlib.pyplot.gcf()
    current_figure.suptitle('Count of publications')
    current_figure.savefig("outputs/pharmgkb/publications.png")

    p = subprocess.Popen("zcat outputs/pharmgkb/Gene_Alleles_Allele.Edge.json.gz | jq .from | sort -u", shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    pharmgkb_genes = [line.strip() for line in p.stdout.readlines()]
    p = subprocess.Popen("zcat outputs/g2p/Gene_Alleles_Allele.Edge.json.gz | jq .from | sort -u", shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    g2p_genes = [line.strip() for line in p.stdout.readlines()]
    sample_df = pandas.DataFrame(upsetplot.from_contents({'pharmgkb': pharmgkb_genes, 'g2p': g2p_genes}))
    upsetplot.plot(sample_df, sort_by="cardinality", sum_over=False, show_counts='%d')
    current_figure = matplotlib.pyplot.gcf()
    current_figure.suptitle('Count of genes')
    current_figure.savefig("outputs/pharmgkb/genes.png")


def transform(find_cmd='ls -1 outputs/pharmgkb/*.*', file_path='/tmp/edges.json', output_dir='outputs/pharmgkb'):
    """Creates schema files."""
    # G = create_graph(load_edges(deduce_edges(find_cmd, file_path)))
    # count_vertexes(find_cmd, G)
    # draw_graph(G, f"{output_dir}/graph.png")
    draw_intersections()


if __name__ == "__main__":
    transform()
