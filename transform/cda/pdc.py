import json
import requests
from attrdict import AttrDict
from string import Template
import yaml
import pathlib
from ioutils import reader



def run_query(query, **kwargs):
    offset = 0
    limit = 1000
    results = {'data': {}}
    while True:
        _query = Template(query).substitute(offset=offset, limit=limit, **kwargs)
        if 'debug' in kwargs:
            print(_query)
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
                    results['data'][entity].append(k)
                else:
                    results['data'][entity].extend(response_entity[k])
        if 'pagination' not in results:
            break
        offset = offset + results['pagination']['count']
    return results


def save_query(obj, query):
    if 'data' in obj:
        entity = list(obj['data'].keys())[0]
        if entity:
            entities = obj['data'][entity]
        else:    
            entities = obj['data']
        with open('source/pdc/{}.json'.format(query), 'w') as outs:
            for e in entities:
                json.dump(e, outs, separators=(',', ':'))
                outs.write('\n')
    else:
        print(query, obj)
        raise Exception('No data')


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
    with open('pdc_queries.yml', 'r') as infile:
        return yaml.load(infile, Loader=yaml.FullLoader)


def paginatedCasesSamplesAliquots(queries=None):
    # expected to run without issue
    if not queries:
        queries = read_queries()
    query_name = 'paginatedCasesSamplesAliquots'
    graphql = queries[query_name]
    save_query(run_query(graphql), query_name)
    

def allCases(queries=None):
    if not queries:
        queries = read_queries()
    query_name = 'allCases'
    graphql = queries[query_name]
    save_query(run_query(graphql), query_name)
      

def allPrograms(queries=None):
    if not queries:
        queries = read_queries()
    query_name = 'allPrograms'
    graphql = queries[query_name]
    save_query(run_query(graphql), query_name)

def getPaginatedFiles(queries=None):
    if not queries:
        queries = read_queries()
    query_name = 'getPaginatedFiles'
    graphql = queries[query_name]
    save_query(run_query(graphql), query_name)

def run_and_save(query_name, queries=None, **kwargs):
    if not queries:
        queries = read_queries()
    graphql = queries[query_name]
    save_query(run_query(graphql, **kwargs), query_name)


def load(name):
    objs = []
    with reader('source/pdc/{}.json'.format(name)) as input:
        for obj in input:
            assert isinstance(obj, dict)
            objs.append(AttrDict(obj))
    assert len(objs) > 1, 'Should have more than 1 record'
    return objs
    
def run_and_load(name, **kwargs):
    run_and_save(name, **kwargs)
    return load(name)



def harvest():
    """Retrieves and saves all relevant data."""
    queries = read_queries()
    pathlib.Path("source/pdc").mkdir(parents=True, exist_ok=True)

    all_cases = run_and_save('allCases')
    print('all_cases', len(all_cases))

    all_programs = run_and_save('allPrograms')
    print('all_programs', len(all_programs))

    all_studies = run_and_load('allStudies')
    print('all_studies', len(all_studies))
    with open('source/pdc/allStudies.json', 'w') as outs:
        for s in all_studies:
            s.biospecimens = run_and_load('biospecimenPerStudy', study_id='"{}"'.format(s.study_id))
            print('biospecimenPerStudy', len(s.biospecimens))            
            s.files = run_and_load('filesPerStudy', study_id='"{}"'.format(s.study_id))
            print('filesPerStudy', len(s.files))            
            json.dump(s, outs, separators=(',', ':'))
            outs.write('\n')

            
def load_all():
    """Returns a tuple with all relevant data."""
    all_cases = load('allCases')
    all_programs = load('allPrograms')
    all_studies = load('allStudies')
    return (all_cases, all_programs, all_studies)        


def old_transform():


    pathlib.Path("source/pdc").mkdir(parents=True, exist_ok=True)

    need_no_parameters = ['studyExperimentalDesign', 'biospecimenPerStudy', 'clinicalPerStudy', 'protocolPerStudy', 'study', 'clinicalMetadata', 'uiSunburstChart', 'quantitiveDataCPTAC2', 'programsProjectsStudies', 'fileMetadata', 'allCases', 'allPrograms', 'tissueSitesAvailable', 'diseasesAvailable', 'allExperimentTypes', 'diseaseTypesPerProject', 'projectsPerExperimentType', 'filesCountPerStudy', 'filesPerStudy', 'projectsPerInstrument', 'uiProtocol', 'pdcDataStats', 'workflowMetadata', 'uiStudy', 'uiCase', 'uiFile', 'uiTissueSiteCaseCount', 'uiAnalyticalFractionsCount', 'uiExperimentBar', 'uiExperimentPie', 'uiPublication']
    needs_parameters = ['experimentalMetadata', 'casePerFile', 'uiExperimentFileCount', 'uiDataCategoryFileCount', 'uiGeneStudySpectralCount', 'uiGeneAliquotSpectralCount', 'getPaginatedFiles']

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
    harvest()
