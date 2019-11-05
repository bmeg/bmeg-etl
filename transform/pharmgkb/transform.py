from bmeg.emitter import JSONEmitter
from fetch_sqlite import fetch_all, sqlite_connection
from bmeg import G2PAssociation, Project, Publication
import json
from bmeg.enrichers.gene_enricher import get_gene
from bmeg.enrichers.phenotype_enricher import phenotype_factory
from bmeg.enrichers.drug_enricher import compound_factory

from bmeg import (G2PAssociation_Publications_Publication, G2PAssociation_Genes_Gene, G2PAssociation_Alleles_Allele,
                  G2PAssociation_Phenotypes_Phenotype, G2PAssociation_Compounds_Compound, G2PAssociation_GenomicFeatures_GenomicFeature,
                  Allele_Gene_Gene)

from bmeg import (Gene, Allele, GenomicFeature)


def create_G2PAssociation(simplified_association, annotation):
    """Create Vertex."""
    association = {}
    association['source_url'] = f"https://www.pharmgkb.org/clinicalAnnotation/{annotation['clinical_annotation_id']}"
    association['evidence_label'] = annotation['level_of_evidence']
    association['description'] = annotation['annotation_text']
    association_parms = {
        'source': 'pharmgkb',
        'source_document': json.dumps(annotation,
                                      sort_keys=True,
                                      separators=(',', ':')),
    }
    for f in ['description', 'evidence_label', 'response_type',
              'oncogenic', 'source_url']:
        association_parms[f] = association.get(f, None)
    return G2PAssociation(
        id=G2PAssociation.make_gid(**association_parms),
        project_id=Project.make_gid("Reference"),
        **association_parms
    )


def create_evidence(hit):
    """Creates Evidence."""
    annotation = hit[hit['source']]
    publications = set([])
    # format evidence as bmeg friendly
    association = hit['association']
    # association['pubmed_ids'] = ','.join([f"https://www.ncbi.nlm.nih.gov/pubmed/{pmid}" ])
    for pmid in annotation['pubmed_ids']:
        url = f"https://www.ncbi.nlm.nih.gov/pubmed/{pmid}"
        publications.add(
            Publication(url=url, title=None, abstract=None, text=None, date=None, author=None, citation=None,
                        id=Publication.make_gid(url), project_id=Project.make_gid("Reference"))
        )
    association.publications = publications
    publication_gids = list(set([p.gid() for p in publications]))
    return (hit, publication_gids)


def drugs(annotation):
    """Selects drug ontology id. (first)"""
    _a = []
    for c in annotation['related_chemicals']:
        for v in c['external_vocabulary']:
            if not v:
                continue
            _a.append(v)
            break
    return _a


def diseases(annotation):
    """Selects disease ontology id. (first)"""
    _a = []
    for c in annotation['related_diseases']:
        for v in c['external_vocabulary']:
            if not v:
                continue
            _a.append(v)
            break
    return _a


PUBLICATIONS = set()
ALLELES = set()
GENOMIC_FEATURES = set()
PHENOTYPES = set()
DRUGS = set()


def toGraph(simplified_association, annotation, emitter):
    """ tuple to graph edges and vertexes """
    association = create_G2PAssociation(simplified_association, annotation)
    association_gid = association.gid()
    emitter.emit_vertex(association)
    # print(association.source_url)

    # assume pubmed transformer creating publication vertex
    publications = []
    for url in simplified_association['publications']:
        publications.append(
            Publication(url=url.strip(), title=None, abstract=None, text=None, date=None, author=None, citation=None,
                        id=Publication.make_gid(url.strip()), project_id=Project.make_gid("Reference"))
        )
    for p in publications:
        publication_gid = p.gid()
        if publication_gid not in PUBLICATIONS:
            emitter.emit_vertex(p)
            PUBLICATIONS.add(publication_gid)
        emitter.emit_edge(
            G2PAssociation_Publications_Publication(
                from_gid=association_gid,
                to_gid=publication_gid
            ),
            emit_backref=True
        )

    # note we assume gene vertexes are already created
    for g in simplified_association['genes']:
        if len(g) > 0:
            emitter.emit_edge(
                G2PAssociation_Genes_Gene(
                    association_gid,
                    Gene.make_gid(g)
                ),
                emit_backref=True
            )

    for feature in simplified_association['features']:
        # print(feature['allele_string'], feature['rsid'])
        allele_strings = feature['allele_string'].split('/')
        ref = allele_strings[0]
        for alt in allele_strings[:1]:
            gene_symbols = feature.get('gene_symbols')
            hugo_symbol = None
            for g in gene_symbols:
                hugo_symbol = g
                break
            allele_gid = Allele.make_gid(
                feature['assembly_name'], feature['seq_region_name'],
                feature['start'], feature['end'],
                ref, alt
            )
            params = {
                'genome': feature['assembly_name'],
                'chromosome': feature['seq_region_name'],
                'start': feature['start'],
                'end': feature['end'],
                'reference_bases': ref,
                'alternate_bases': alt,
                'strand': '+',
                'hugo_symbol': hugo_symbol,
                'id': allele_gid,
                'project_id': Project.make_gid("Reference")
            }
            allele = Allele(**params)
            allele.validate()
            if allele_gid not in ALLELES:
                ALLELES.add(allele_gid)
                emitter.emit_vertex(allele)
                if hugo_symbol:
                    try:
                        g = get_gene(hugo_symbol)
                        gene_gid = Gene.make_gid(g['ensembl_gene_id'])
                        emitter.emit_edge(
                            Allele_Gene_Gene(
                                allele_gid,
                                gene_gid,
                            ),
                            emit_backref=True
                        )
                    except Exception as e:
                        print('WARNING Gene.make_gid', e, hugo_symbol, g)

            emitter.emit_edge(
                G2PAssociation_Alleles_Allele(
                    association_gid,
                    allele_gid
                ),
                emit_backref=True
            )

    for feature_name in simplified_association['feature_names']:
        gf_id = GenomicFeature.make_gid(
            None, None,
            None, None,
            None, feature_name
        )
        params = {
            'name': feature_name,
            'id': gf_id,
            'project_id': Project.make_gid("Reference")
        }
        gf = GenomicFeature(**params)
        gf.validate()
        if gf_id not in GENOMIC_FEATURES:
            emitter.emit_vertex(gf)
            GENOMIC_FEATURES.add(gf_id)
        emitter.emit_edge(
            G2PAssociation_GenomicFeatures_GenomicFeature(
                association_gid,
                gf_id
            ),
            emit_backref=True
        )

    for p in simplified_association['diseases']:
        p = phenotype_factory(p['name'])
        p_gid = p.gid()
        if p_gid not in PHENOTYPES:
            emitter.emit_vertex(p)
            PHENOTYPES.add(p_gid)
        emitter.emit_edge(
            G2PAssociation_Phenotypes_Phenotype(
                association_gid,
                p_gid
            ),
            emit_backref=True
        )

    for d in simplified_association['drugs']:
        d = compound_factory(d['name'])
        d_gid = d.gid()
        if d_gid not in DRUGS:
            emitter.emit_vertex(d)
            DRUGS.add(d_gid)
        emitter.emit_edge(
            G2PAssociation_Compounds_Compound(
                association_gid,
                d_gid
            ),
            emit_backref=True
        )


def transform(source_path="source/pharmgkb",
              emitter_prefix=None,
              emitter_directory='pharmgkb'):
    """"""
    emitter = JSONEmitter(directory=emitter_directory, prefix=emitter_prefix)
    # clinical_ann_metadata
    conn = sqlite_connection(source_path)
    for association in fetch_all(conn=conn):
        # fix genes as sometimes they have a CSV field
        genes = []
        for ensembl_ids in [g['ensembl_id'] for g in association['gene']]:
            for ensembl_id in ensembl_ids.split(','):
                genes.append(ensembl_id.replace('"', '').strip())
        simplified_association = {
            'level_of_evidence': association['level_of_evidence'],
            'genes': genes,
            'diseases': diseases(association),
            'drugs': drugs(association),
            'feature_names': association['feature_names'],
            'features': association['features'],
            'description': ' '.join([v['sentence'] for v in association['variant_annotations_ids']]),
            'publications': [f"https://www.ncbi.nlm.nih.gov/pubmed/{pmid}" for pmid in association['pubmed_ids']]
        }
        toGraph(simplified_association, association, emitter)
        # if association['clinical_annotation_id'] == '1447989723':
        #     print(json.dumps(simplified_association, separators=(',',':')))
    emitter.close()


if __name__ == '__main__':  # pragma: no cover
    transform()
