
import json
from bmeg.vertex import G2PAssociation


def create_association(hit):
    """ create Vertex """
    association = hit['association']
    association_parms = {
        'source': hit['source'],
        'source_document': json.dumps(hit[hit['source']],
                                      separators=(',', ':')),
    }
    for f in ['description', 'evidence_label', 'response_type',
              'oncogenic', 'source_url']:
        association_parms[f] = association.get(f, None)
    return G2PAssociation(**association_parms)


def normalize(hit):
    """ returns a tuple of (hit, association), where hit association
    has been modified to remove edge keys and gid inserted
    """
    association = create_association(hit)
    return (hit, association)
