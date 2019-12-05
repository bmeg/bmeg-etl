

import json
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
from bmeg import (G2PAssociation_Publications_Publication, G2PAssociation_Genes_Gene, G2PAssociation_Alleles_Allele,
                  G2PAssociation_Phenotypes_Phenotype, G2PAssociation_Compounds_Compound, G2PAssociation_GenomicFeatures_GenomicFeature,
                  GenomicFeature_Genes_Gene)
from bmeg import Deadletter
from bmeg.emitter import new_emitter


def normalizeAssociations(path):
    """ create a record from input """
    input_stream = bmeg.ioutils.reader(path)
    # create a tuple where:
    # vertices - dict of vertices to be created
    # genes ... publications - arrays of edges from association
    NormalizedAssociation = collections.namedtuple(
        'NormalizedAssociation',
        ['vertices', 'genes', 'features', 'environments', 'phenotypes',
         'publications', 'association', 'genomic_features', 'missing_vertexes']
    )
    for line in input_stream:
        hit = json.loads(line)
        if hit['source'] == 'litvar':
            continue
        (hit, genes, missing_genes) = genes_normalize(hit)
        (hit, features, genomic_features, missing_features) = features_normalize(hit)
        (hit, environments) = environments_normalize(hit)
        (hit, phenotypes) = phenotypes_normalize(hit)
        (hit, publications) = publication_normalize(hit)
        (hit, association) = association_normalize(hit)
        yield NormalizedAssociation(hit, genes, features, environments,
                                    phenotypes, publications, association,
                                    genomic_features,
                                    missing_genes + missing_features)


def toGraph(normalized_association, emitter):
    """ tuple to graph edges and vertexes """
    na = normalized_association

    association = na.association
    emitter.emit_vertex(association)

    # assume pubmed transformer creating publication vertex
    for publication_gid in na.publications:
        emitter.emit_edge(
            G2PAssociation_Publications_Publication(
                from_gid=association.gid(),
                to_gid=publication_gid
            ),
            emit_backref=True
        )

    # note we assume gene vertexes are already created
    for gene_gid in na.genes:
        emitter.emit_edge(
            G2PAssociation_Genes_Gene(
                association.gid(),
                gene_gid
            ),
            emit_backref=True
        )

    for allele in na.vertices['features']:
        emitter.emit_vertex(allele)
    for feature_gid in na.features:
        emitter.emit_edge(
            G2PAssociation_Alleles_Allele(
                association.gid(),
                feature_gid
            ),
            emit_backref=True
        )

    for feature in na.vertices['genomic_features']:
        emitter.emit_vertex(feature)
    for feature_gid in na.genomic_features:
        emitter.emit_edge(
            G2PAssociation_GenomicFeatures_GenomicFeature(
                association.gid(),
                feature_gid
            ),
            emit_backref=True
        )

    for genomic_feature_has_gene in na.vertices['genomic_feature_has_gene']:
        emitter.emit_edge(
            GenomicFeature_Genes_Gene(
                genomic_feature_has_gene[0],
                genomic_feature_has_gene[1],
            ),
            emit_backref=True
        )

    for phenotype in na.vertices['phenotypes']:
        emitter.emit_vertex(phenotype)
    for phenotype_gid in na.phenotypes:
        emitter.emit_edge(
            G2PAssociation_Phenotypes_Phenotype(
                association.gid(),
                phenotype_gid
            ),
            emit_backref=True
        )

    for environment in na.vertices['environments']:
        emitter.emit_vertex(environment)
    for environment_gid in na.environments:
        emitter.emit_edge(
            G2PAssociation_Compounds_Compound(
                association.gid(),
                environment_gid
            ),
            emit_backref=True
        )

    EMITTED_DEADLETTER = {}
    for m in na.missing_vertexes:
        dl = Deadletter(**m)
        if dl.gid() in EMITTED_DEADLETTER:
            continue
        emitter.emit_vertex(dl)
        EMITTED_DEADLETTER[dl.gid()] = True


def transform(input_path, emitter_directory, emitter_class):
    """ parse the association and write to graph using emitter"""
    emitter = new_emitter(name=emitter_class, directory=emitter_directory)
    association_cache = []
    for normalized_association in normalizeAssociations(input_path):
        if normalized_association.association.gid() in association_cache:
            continue
        association_cache.append(normalized_association.association.gid())
        toGraph(normalized_association, emitter)
    emitter.close()


def main():  # pragma: no cover
    parser = default_argument_parser(emitter_directory_default='g2p')
    parser.add_argument('--input_path', type=str,
                        default='source/g2p/all.json.gz',
                        help='Path to g2p associations for import')
    options = parser.parse_args()
    default_logging(options.loglevel)
    transform(options.input_path,
              emitter_directory=options.emitter_directory,
              emitter_class=options.emitter)


if __name__ == '__main__':
    main()
