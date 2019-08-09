import json
import requests
from attrdict import AttrDict
from string import Template
import yaml
import pathlib
from transform.pdc.ioutils import reader


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
                if isinstance(k, dict):
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


def read_queries():
    with open('transform/pdc/pdc_queries.yml', 'r') as infile:
        return yaml.load(infile, Loader=yaml.FullLoader)


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
    read_queries()
    pathlib.Path("source/pdc").mkdir(parents=True, exist_ok=True)

    run_and_save('allCases')

    run_and_save('allPrograms')

    all_studies = run_and_load('allStudies')
    with open('source/pdc/allStudies.json', 'w') as outs:
        for s in all_studies:
            s.biospecimens = run_and_load('biospecimenPerStudy', study_id='"{}"'.format(s.study_id))
            s.files = run_and_load('filesPerStudy', study_id='"{}"'.format(s.study_id))
            json.dump(s, outs, separators=(',', ':'))
            outs.write('\n')


def load_all():
    """Returns a tuple with all relevant data."""
    all_cases = load('allCases')
    all_programs = load('allPrograms')
    all_studies = load('allStudies')
    return (all_cases, all_programs, all_studies)


if __name__ == "__main__":
    harvest()
