import json
import requests
from attrdict import AttrDict
from string import Template
import yaml
import pathlib


# see https://pdc.esacinc.com/data-dictionary/publicapi-documentation/#!/Paginated_Records/getPaginatedCases_0
# https://pdc.esacinc.com/data-dictionary/harmonization.html

def get_pdc_schema():
    """Returns (entities, queries)"""

    introspection_query = """

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
    r = requests.get('https://pdc.esacinc.com/graphql?query={}'.format(introspection_query))
    pdc = AttrDict(r.json())

    # xform
    entities = {t['name']: [f for f in t['fields']] for t in pdc.data['__schema']['types'] if t['kind'] == 'OBJECT' and not t['name'].startswith('_')}
    queries = {q['name']:q  for q in next(filter(lambda t: t['name'] == 'Query', pdc.data['__schema']['types']))['fields']}
    return (entities, queries)
    # print(entities.keys())
    # print(queries.keys())





def sliceindex(x):
    i = 0
    for c in x:
        if c.isalpha():
            i = i + 1
            return i
        i = i + 1

def upperfirst(x):
    i = sliceindex(x)
    return x[:i].upper() + x[i:]


def make_query(query, entity):
    fields = [f['name'] for f in [t['fields'] for t in pdc.data['__schema']['types'] if t['name'] == entity][0]]
    q = '{{ {} {{ {} }} }}'.format(query, ' '.join(sorted(fields)))
    print(q)
    return q


def run_query(query):
    offset = 0
    limit = 1000
    results = {'data': {}}
    while True:
        _query = Template(query).substitute(offset=offset, limit=limit)
        r = requests.get('https://pdc.esacinc.com/graphql?query={}'.format(_query))
        response = r.json()
        if 'errors' in response:
            return response
        if 'data' in response:
            entity = list(response['data'].keys())[0]
            if entity not in results['data']:
               results['data'][entity] = []
            response_entity = response['data'][entity]
            for k in response_entity:
                if k in ['total', 'pagination']:
                    continue
                if isinstance(k,dict):
                    results['data'][entity].extend(k)
                else:
                    results['data'][entity].extend(response_entity[k])
        if 'pagination' not in results:
            break
        offset = offset + results['pagination']['count']
    return results


def save_query(obj, query):
    if 'data' in obj:
        entity = list(obj['data'].keys())[0]
        print(entity, len(obj['data'][entity]))
    else:
        print(query, obj)
        raise Exception('No data')
    with open('source/pdc/{}.json'.format(query), 'w') as outs:
        json.dump(obj, outs)


def parse_query_name(query):
    parts = query.split(' ')
    if len(parts) < 2:
        return None
    return parts[1]


def query_return_types(entities, queries):
    # get expected return type and it's attributes
    for k, v in queries.items():
        v = AttrDict(v)
        if 'ofType' in v.type and v.type.ofType:
            print(v.type.ofType.name, len(sorted(entities[v.type.ofType.name])) )

def field_format(field,type, entities):
    if type['kind'] == 'OBJECT':
        typed_fields = {e['name']:e for e in entities[type['name']] if e['type']['kind'] == 'SCALAR'}
        return '\n    {} {{ {} }}\n'.format(field,field_list(typed_fields, entities))
    if type['kind'] == 'LIST':
        list_type = type['ofType']['name']
        typed_fields = sorted([e['name'] for e in entities[list_type] if e['type']['kind'] == 'SCALAR'])
        return '\n    {} {{ {} }}\n'.format(field, ' '.join(typed_fields))
    return '\n    {}'.format(field)


def field_list(typed_fields, entities):
    fields = sorted([e for e in typed_fields])
    types = [typed_fields[f]['type'] for f in fields]
    return ' '.join([field_format(field,type, entities) for field,type in zip(fields, types)])


def naive_queries(entities=None, queries=None):
    """Generates dumb queries"""
    if not entities:
        (entities, queries) = get_pdc_schema()
    query_strings = {}
    for k, v in queries.items():
        v = AttrDict(v)
        if 'ofType' in v.type and v.type.ofType:
            typed_fields = {e['name']:e for e in entities[v.type.ofType.name]}
            q = '    {{\n    {} {{ {}\n    }} \n    }}'.format(k, field_list(typed_fields, entities))
            query_strings[k] = q
    return AttrDict(query_strings)


# create basic queries, write queries to editiable yaml file

def generate_queries():
    naive_queries = pdc.naive_queries()
    with open('pdc_queries.yml', 'w') as outfile:
        for k, v in naive_queries.items():
            outfile.write('{}: |\n'.format(k))
            outfile.write('{}\n\n'.format(v))

# read yaml files

def read_queries():
    with open('transform/pdc/pdc_queries.yml', 'r') as infile:
        return yaml.load(infile, Loader=yaml.FullLoader)


def transform():


    pathlib.Path("source/pdc").mkdir(parents=True, exist_ok=True)

    need_no_parameters = ['studyExperimentalDesign', 'biospecimenPerStudy', 'clinicalPerStudy', 'protocolPerStudy', 'study', 'clinicalMetadata', 'uiSunburstChart', 'quantitiveDataCPTAC2', 'programsProjectsStudies', 'fileMetadata', 'allCases', 'allPrograms', 'tissueSitesAvailable', 'diseasesAvailable', 'allExperimentTypes', 'diseaseTypesPerProject', 'projectsPerExperimentType', 'filesCountPerStudy', 'filesPerStudy', 'projectsPerInstrument', 'uiProtocol', 'pdcDataStats', 'workflowMetadata', 'uiStudy', 'uiCase', 'uiFile', 'uiTissueSiteCaseCount', 'uiAnalyticalFractionsCount', 'uiExperimentBar', 'uiExperimentPie', 'uiPublication']
    needs_parameters = ['experimentalMetadata', 'casePerFile', 'uiExperimentFileCount', 'uiDataCategoryFileCount', 'uiGeneStudySpectralCount', 'uiGeneAliquotSpectralCount']

    queries = read_queries()

    # expected to run without issue
    for query_name in need_no_parameters:
        graphql = queries[query_name]
        try:
            save_query(run_query(graphql), query_name)
        except Exception as e:
            print(query_name, e)

    # expected to run without issue
    query_name = 'paginatedCasesSamplesAliquots'
    graphql = queries[query_name]
    save_query(run_query(graphql), query_name)

    print('>>>> expected to have errors <<<<')
    for query_name in needs_parameters:
        graphql = queries[query_name]
        try:
            save_query(run_query(graphql), query_name)
        except Exception as e:
            print(query_name, e)

if __name__ == "__main__":
    transform()
