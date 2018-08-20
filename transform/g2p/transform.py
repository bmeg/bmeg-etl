

import json
import sys
import collections

from transform.g2p.genes import normalize as genes_normalize
from transform.g2p.features import normalize as features_normalize
from transform.g2p.environments import normalize as environments_normalize
from transform.g2p.phenotypes import normalize as phenotypes_normalize
from transform.g2p.association import normalize as association_normalize
from transform.g2p.publications import normalize as publication_normalize

import bmeg.ioutils
from bmeg.util.logging import default_logging
from bmeg.util.cli import default_argument_parser

from bmeg.emitter import JSONEmitter as Emitter

files = {}


def normalizeAssociations(path):
    """ create a record from input """
    input_stream = bmeg.ioutils.reader(path)
    NormalizedAssociation = collections.namedtuple(
        'NormalizedAssociation',
        ['vertices', 'genes', 'features', 'environments', 'phenotypes',
         'publications', 'association'])
    for line in input_stream:
        hit = json.loads(line)
        (hit, genes) = genes_normalize(hit)
        (hit, features) = features_normalize(hit)
        (hit, environments) = environments_normalize(hit)
        (hit, phenotypes) = phenotypes_normalize(hit)
        (hit, publications) = publication_normalize(hit)
        (hit, association) = association_normalize(hit)
        yield NormalizedAssociation(hit, genes, features, environments,
                                    phenotypes, publications, association)


def writeGraph(cls, obj):
    """ write Vertex and Edges """
    pass
    # TODO: emitter ...


def toGraph(normalized_association, prefix):
    """ tuple to graph edges and vertexes """
    na = normalized_association
    association = na.association
    for gene_gid in na.vertices['genes']:
        print('create edge {}->{}'.format(association.gid(), gene_gid))
    for feature_gid in na.vertices['features']:
        print('create edge {}->{}'.format(association.gid(), feature_gid))
    for environment_gid in na.environments:
        print('create edge {}->{}'.format(association.gid(), environment_gid))
    for phenotype_gid in na.phenotypes:
        print('create edge {}->{}'.format(association.gid(), phenotype_gid))
    for publication_gid in na.publications:
        print('create edge {}->{}'.format(association.gid(), publication_gid))
    print('create vertex {}'.format(association.gid()))

    # association_gid = normalized_association.association['gid']

    # # write genes
    # for gid in normalized_association.genes.keys():
    #     # write data
    #     data = dict(genes[gid])
    #     data['id'] = gid
    #     writeGraph('Gene', data)
    #
    # # write features
    # for gid in normalized_association.features.keys():
    #     #  write data
    #     data = features[gid]
    #     if not gid.startswith('gene'):
    #         data['id'] = gid
    #         writeGraph('Variant', data)
    #
    # # write environments
    # for gid in normalized_association.environments.keys():
    #     data = environments[gid]
    #     data['id'] = gid
    #     writeGraph('Compound', data)
    #
    # # write phenotypes
    # for gid in normalized_association.phenotypes.keys():
    #     data = phenotypes[gid]
    #     data['id'] = gid
    #     writeGraph('Phenotype', data)
    #
    # # write association data
    # del association['gid']
    # association['id'] = association_gid
    # writeGraph('G2PAssociation', association)


def transform(input_path, prefix):
    """ parse the association and write to graph """
    for normalized_association in normalizeAssociations(input_path):
        toGraph(normalized_association, prefix)


def main():  # pragma: no cover
    parser = default_argument_parser()
    parser.add_argument('--input_path', type=str,
                        default='source/g2p/all.json.gz',
                        help='Path to g2p associations for import')

    # We don't need the first argument, which is the program name
    options = parser.parse_args(sys.argv[1:])
    default_logging(options.loglevel)
    transform(options)


if __name__ == '__main__':
    main()
