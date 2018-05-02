import sys
import json
import argparse

import convert.sample_pb2 as schema
import convert.record as record

# fetch data ------------
# python gdc-scan.py cases list

def sample_gid(data):
    return 'biosample:' + data

class CohortGenerator(record.RecordGenerator):
    def __init__(self):
        super(CohortGenerator, self).__init__('Cohort')

    def schema(self):
        return schema.Cohort()

    def gid(self, data):
        return 'cohort:' + data['project']['project_id']

    def update(self, cohort, data):
        if not cohort.name:
            cohort.name = data['project']['project_id']
        if not cohort.location:
            cohort.location = data['project']['primary_site']
        if not cohort.description:
            cohort.description = data['project']['name']

        if 'submitter_sample_ids' in data:
            for sample in data['submitter_sample_ids']:
                record.append_unique(cohort.hasMemberEdges, sample_gid(sample))

        return cohort

def read_input(path):
    with open(path) as input:
        return input.read().split('\n')

def convert_cohort(input, output):
    state = {
        'Cohort': {},
        'types': ['Cohort'],
        'generators': {
            'cohort': CohortGenerator()
        }
    }

    lines = read_input(input)
    for line in lines:
        if not line == '':
            try:
                sample = json.loads(line)
            except:
                print('sample failed! ' + line)

            state['generators']['cohort'].find(state, sample)

    record.output_state(state, output)

def parse_args(args):
    args = args[1:]
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--input', type=str, help='path to input file')
    parser.add_argument('--output', type=str, help='path to output file')
    
    return parser.parse_args(args)

if __name__ == '__main__':
    options = parse_args(sys.argv)
    convert_cohort(options.input, options.output)