import json
import requests
from attrdict import AttrDict
from string import Template
import yaml
import pathlib
from ioutils import reader

# see https://pdc.esacinc.com/data-dictionary/publicapi-documentation/#!/Paginated_Records/getPaginatedCases_0
# https://pdc.esacinc.com/data-dictionary/harmonization.html

INTROSPECTION_QUERY = """
    query IntrospectionQuery {
      __schema {
        queryType { name }
        mutationType { name }
        subscriptionType { name }
        types {
          ...FullType
        }
        directives {
          name
          description
          locations
          args {
            ...InputValue
          }
        }
      }
    }

    fragment FullType on __Type {
      kind
      name
      description
      fields(includeDeprecated: true) {
        name
        description
        args {
          ...InputValue
        }
        type {
          ...TypeRef
        }
        isDeprecated
        deprecationReason
      }
      inputFields {
        ...InputValue
      }
      interfaces {
        ...TypeRef
      }
      enumValues(includeDeprecated: true) {
        name
        description
        isDeprecated
        deprecationReason
      }
      possibleTypes {
        ...TypeRef
      }
    }

    fragment InputValue on __InputValue {
      name
      description
      type { ...TypeRef }
      defaultValue
    }

    fragment TypeRef on __Type {
      kind
      name
      ofType {
        kind
        name
        ofType {
          kind
          name
          ofType {
            kind
            name
            ofType {
              kind
              name
              ofType {
                kind
                name
                ofType {
                  kind
                  name
                  ofType {
                    kind
                    name
                  }
                }
              }
            }
          }
        }
      }
    }

"""

TEMPLATE = """{{
  {} ({})
  {{
    {}
  }}
}}
"""

class PDCSchema():
    """Introspects PDC schema, provides useful functions."""    

    def __init__(self):
        r = requests.get('https://pdc.esacinc.com/graphql?query={}'.format(INTROSPECTION_QUERY))    
        self.schema = AttrDict(r.json())

    def queries_dict(self):
        """Dict of queries"""
        query_list = [t['fields'] for t in self.schema.data['__schema']['types'] if t['name'] == 'Query'][0]
        queries = {q['name']:q for q in query_list}
        return queries

    def query_fields(self, name):
        """Fields returned from query"""
        return self.entity_fields(self.queries_dict()[name]['type']['ofType']['name'])

    def query_args(self, name):
        """Arguments to a query"""
        query_list = [t['fields'] for t in self.schema.data['__schema']['types'] if t['name'] == 'Query'][0]
        queries = {q['name']:q for q in query_list}
        return sorted([a['name'] for a in self.queries_dict()[name]['args']])

    def entity_fields(self, name):
        """Fields returned from entity"""
        o = [t for t in self.schema.data['__schema']['types'] if t['name'] == name][0]
        return sorted([f['name'] for f in o['fields']])

    def entities_dict(self):
        """Dict of entities"""
        entity_list = [t for t in self.schema.data['__schema']['types']]
        entities = {e['name']:e for e in entity_list}
        return entities

    def make_query(self, name):
        """Naive query"""
        return TEMPLATE.format(name,
                        ' '.join(['{}:${}'.format(a,a) for a in self.query_args(name)]),
                        '\n    '.join(self.query_fields(name)))    
        
