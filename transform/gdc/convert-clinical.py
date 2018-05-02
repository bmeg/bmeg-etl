#!/usr/bin/env python

# download:
# ../gdc_scan.py files download --project TCGA-LUAD --type "Clinical Supplement"
# example usage:
# python -m convert.gdc.convert-clinical Individual --input ~/Data/gdc/clinical/ --output ~/Data/gdc/individual.json

import os
import re
import json
import tempfile
import subprocess
from copy import deepcopy
import argparse
from xml.dom.minidom import parseString
from google.protobuf import json_format

#from ga4gh import bio_metadata_pb2
# from clinical_pb2 import Individual, Cohort, Biosample, IndividualCohort, DrugTherapy, RadiationTherapy, ClinicalFollowup

from bmeg import conversions
from bmeg.clinical_pb2 import Individual, Cohort, Biosample, IndividualCohort, DrugTherapy, RadiationTherapy, ClinicalFollowup

def getText(nodelist):
    rc = []
    for node in nodelist:
        if node.nodeType == node.TEXT_NODE:
            rc.append(node.data)
    return ''.join(rc)

def dom_scan(node, query):
    stack = query.split('/')
    if node.localName == stack[0]:
        return dom_scan_iter(node, stack[1:], [stack[0]])

def dom_scan_iter(node, stack, prefix):
    if len(stack):
        for child in node.childNodes:
                if child.nodeType == child.ELEMENT_NODE:
                    if child.localName == stack[0]:
                        for out in dom_scan_iter(child, stack[1:], prefix + [stack[0]]):
                            yield out
                    elif '*' == stack[0]:
                        for out in dom_scan_iter(child, stack[1:], prefix + [child.localName]):
                            yield out
    else:
        if node.nodeType == node.ELEMENT_NODE:
            yield node, prefix, dict(node.attributes.items()), getText( node.childNodes )
        elif node.nodeType == node.TEXT_NODE:
            yield node, prefix, None, getText( node.childNodes )

def force_parse(s):
    try:
        return int(s)
    except ValueError:
        try:
            return float(s)
        except ValueError:
            return s

def record_initial_state(generators):
    state = {}
    state['types'] = generators.keys()
    state['generators'] = generators
    for key in generators:
        state[key] = {}

    return state

def process_input(state, path, source, process):
    if os.path.isdir(path):
        for item in os.listdir(path):
            elem = os.path.join(path, item)
            if os.path.isdir(elem):
                process_input(state, elem, source, process)
            else:
                with open(os.path.join(path, item)) as input:
                    state = process(state, source, input)
    else:
        with open(path) as input:
            state = process(state, source, input)

    return state

def message_to_json(message):
    msg = json_format.MessageToDict(message)
    return json.dumps(msg)

def resolve_type(existing, novel):
    if existing is None:
        return novel
    elif existing == int:
        return novel
    elif existing == float:
        if novel == str or novel == unicode:
            return unicode
        else:
            return float
    else:
        return unicode

def align_info(keys, emit):
    for key in keys:
        into = keys[key]
        if key in emit.attributes and len(emit.attributes[key]) > 0 and not into == str and not into == unicode:
            pass
            #translated = map(into, emit.attributes.attr[key].values)
            #del emit.attributes.attr[key]
            #emit.attributes[key].extend(translated)

    return emit

def output_state(state, path):
    for type in state['types']:
        json = []
        # if type == 'Individual':
        #     print(state['generators']['Individual'].infokeys)
        for key in state[type]:
            emit = state[type][key]
            #if type == 'Individual':
            #    emit = align_info(state['generators']['Individual'].infokeys, emit)
            message = message_to_json(emit)
            json.append(message)
        out = '\n'.join(json)
        # if len(json):

        outpath = 'tcga.' + type + '.json'
        outhandle = open(os.path.join(path, outpath), 'w')
        outhandle.write(out)
        outhandle.close()

class RecordGenerator(object):
    def __init__(self, name):
        self.name = name

    def schema(self):
        raise Exception('schema() not implemented')

    def gid(self, data):
        raise Exception('gid() not implemented')

    def update(self, record, data):
        raise Exception('update() not implemented')

    def create(self, data):
        record = self.schema()
        gid = self.gid(data)
        record.id = gid
        #record.type = self.name
        self.update(record, data)
        return record

    def find(self, data):
        gid = self.gid(data)
        record = state[self.name].get(gid)
        if record is None:
            record = self.create(data)
            state[self.name][gid] = record
        else:
            self.update(record, data)
        return record

class ProjectGenerator(RecordGenerator):
    def __init__(self, individual_gid):
        super(ProjectGenerator, self).__init__('Project')
        self.individual_gid = individual_gid

    def schema(self):
        return IndividualCohort() #matrix_pb2.Project()

    def gid(self, data):
        dataset = data['project_code'] + '-' + data['disease_code']
        return 'individualCohort:' + dataset

    def update(self, cohort, data):
        cohort.name = data['project_code'] + '-' + data['disease_code']
        return cohort

class IndividualGenerator(RecordGenerator):
    def __init__(self):
        super(IndividualGenerator, self).__init__('Individual')
        self.infokeys = {}

    def schema(self):
        return Individual()

    def gid(self, data):
        #dataset = data['project_code'] + '-' + data['disease_code']
        #return 'individual:' + dataset + ':' + data['bcr_patient_barcode']
        return data['bcr_patient_uuid'].lower()

    def update(self, individual, data):
        dataset = data['project_code'] + '-' + data['disease_code']
        individual.source = data['source']
        individual.name = data['bcr_patient_barcode']
        individual.dataset_id = dataset
        if 'submitted_tumor_site' in data:
            individual.attributes['tumorSite'] = data['submitted_tumor_site']

        for key in data:
            if len(data[key]):
                keytype = type(force_parse(data[key]))
                self.infokeys[key] = resolve_type(self.infokeys.get(key), keytype)
                individual.attributes[key] = data[key]

        return individual

class BiosampleGenerator(RecordGenerator):
    def __init__(self, individual_gid):
        super(BiosampleGenerator, self).__init__('Biosample')
        self.individual_gid = individual_gid

    def schema(self):
        return Biosample()

    def gid(self, data):
        #dataset = data['project_code'] + '-' + data['disease_code']
        #return 'biosample:' + dataset + ':' + data['bcr_sample_barcode']
        return data['bcr_sample_uuid'].lower()

    def update(self, sample, data):
        dataset = data['project_code'] + '-' + data['disease_code']
        sample.source = data['source']
        sample.id = data['bcr_sample_uuid'].lower()
        sample.name = data['bcr_sample_barcode']
        sample.dataset_id = dataset
        if 'submitted_tumor_site' in data:
            site = data['submitted_tumor_site']
            sample.disease.term = site

        # sample.sampleType = 'tumor' if re.search('tumor', data['sample_type'], re.IGNORECASE) else 'normal'
        for key in data:
            if len(data[key]):
                sample.attributes[key] = data[key]

        individual_id = {
            'project_code': data['project_code'],
            'disease_code': data['disease_code'],
            'bcr_patient_barcode': data['bcr_sample_barcode'][:12],
            'bcr_patient_uuid' : data['bcr_patient_uuid']
        }

        sample.individual_id = self.individual_gid(individual_id)
        if sample.individual_id in state['Individual']:
            individual = state['Individual'][sample.individual_id]
            if 'tumorSite' in individual.info:
                if len(individual.info['tumorSite']) > 0:
                    site = individual.info['tumorSite'][0]
                    sample.disease.term = site.lower()
        # record.append_unique(sample.sampleOf, self.individual_gid({
        #     'bcr_patient_barcode': data['bcr_sample_barcode'][:12]}))

        return sample

class DrugGenerator(RecordGenerator):
    def __init__(self, id):
        super(DrugGenerator, self).__init__('DrugTherapy')
        self.id = id
        self.pubchem = {}

    def gid(self, data):
        #print data
        #dataset = data['project_code'] + '-' + data['disease_code']
        # return 'drug_therapy:' + data['bcr_drug_barcode']
        return data['bcr_drug_barcode']

    def schema(self):
        return DrugTherapy() #matrix_pb2.Project()

    def update(self, drug, data):
        print data
        drug.source = data['source']
        drug.id = data['bcr_drug_uuid'].lower()
        drug.name = data['bcr_drug_barcode']
        drug.individual_id = self.id(data)
        drug.drug_name = data['pharmaceutical_therapy_drug_name']
        drug.pubchem_id = conversions.pubchem(drug.drug_name.lower())
        # if drug.drug_name.lower().strip() in self.pubchem:
        #     drug.pubchem_id = self.pubchem[drug.drug_name.lower()]
        drug.prescribed_dose = data['prescribed_dose']
        return drug

class RadiationGenerator(RecordGenerator):
    def __init__(self, id):
        super(RadiationGenerator, self).__init__('RadiationTherapy')
        self.id = id

    def gid(self, data):
        print data
        #dataset = data['project_code'] + '-' + data['disease_code']
        # return 'radiation_therapy:' + data['bcr_radiation_barcode']
        return data['bcr_radiation_barcode']

    def schema(self):
        return RadiationTherapy() #matrix_pb2.Project()

    def update(self, rad, data):
        rad.source = data['source']
        rad.id = data['bcr_radiation_uuid'].lower()
        rad.name = data['bcr_radiation_barcode']
        rad.individual_id = self.id(data)
        try:
            v = float(data['radiation_total_dose'])
            rad.total_dose = v
        except ValueError:
            pass
        return rad

class FollowupGenerator(RecordGenerator):
    def __init__(self, id):
        super(FollowupGenerator, self).__init__('ClinicalFollowup')
        self.id = id

    def gid(self, data):
        #print self.id, data
        #dataset = data['project_code'] + '-' + data['disease_code']
        # return 'clinical_followup:' + data['bcr_followup_barcode']
        return data['bcr_followup_barcode']

    def schema(self):
        return ClinicalFollowup() #matrix_pb2.Project()

    def update(self, followup, data):
        followup.source = data['source']
        followup.id = data['bcr_followup_uuid'].lower()
        followup.name = data['bcr_followup_barcode']
        followup.individual_id = self.id(data)
        followup.date = "%s-%s-%s" % (data['form_completion_year'], data['form_completion_month'], data['form_completion_day'])
        followup.vital_status = data['vital_status']
        try:
            days_to_death = int(data['death_days_to'])
            followup.days_to_death = days_to_death
        except ValueError:
            pass
        return followup


def extract_attribute(data, stack, attr, text):
    if 'xsd_ver' in attr:
        p_name = attr.get('preferred_name', stack[-1])
        if len(p_name) == 0:
            p_name = stack[-1]
        data[p_name] = text.strip()

class ClinicalParser:
    def __init__(self, pubchem={}):
        self.pubchem = pubchem

    def parseXMLFile(self, state, source, dom):
        root_node = dom.childNodes[0]
        admin = {}
        for node, stack, attr, text in dom_scan(root_node, 'tcga_bcr/admin/*'):
            admin[stack[-1]] = text.strip()

        patient_barcode = None
        for node, stack, attr, text in dom_scan(root_node, 'tcga_bcr/patient/bcr_patient_barcode'):
            patient_barcode = text.strip()

        patient_data = admin
        patient_data['source'] = source
        for node, stack, attr, text in dom_scan(root_node, 'tcga_bcr/patient/*'):
            extract_attribute(patient_data, stack, attr, text)

        found = False
        for node, stack, attr, text in dom_scan(root_node, 'tcga_bcr/patient/stage_event/*'):
            extract_attribute(patient_data, stack, attr, text)
            found = True
        for node, stack, attr, text in dom_scan(root_node, 'tcga_bcr/patient/stage_event/*/*'):
            extract_attribute(patient_data, stack, attr, text)
            found = True
        for node, stack, attr, text in dom_scan(root_node, 'tcga_bcr/patient/stage_event/tnm_categories/*/*'):
            extract_attribute(patient_data, stack, attr, text)
            found = True
        if found:
            self.emit( state, patient_barcode, patient_data, "Individual" )

        for s_node, s_stack, s_attr, s_text in dom_scan(root_node, 'tcga_bcr/patient/samples/sample'):
            sample_barcode = None
            for c_node, c_stack, c_attr, c_text in dom_scan(s_node, 'sample/bcr_sample_barcode'):
                sample_barcode = c_text.strip()
            sample_data = deepcopy(patient_data)
            sample_data['source'] = source
            for c_node, c_stack, c_attr, c_text in dom_scan(s_node, 'sample/*'):
                if 'xsd_ver' in c_attr:
                    sample_data[c_attr.get('preferred_name', c_stack[-1])] = c_text.strip()
            self.emit( state, sample_barcode, sample_data, "Biosample" )

        """
        for s_node, s_stack, s_attr, s_text in dom_scan(root_node, 'tcga_bcr/patient/samples/sample/portions/portion'):
            portion_barcode = None
            for c_node, c_stack, c_attr, c_text in dom_scan(s_node, 'portion/bcr_portion_barcode'):
                portion_barcode = c_text
            portion_data = {}
            for c_node, c_stack, c_attr, c_text in dom_scan(s_node, 'portion/*'):
                if 'xsd_ver' in c_attr:
                    portion_data[c_attr.get('preferred_name', c_stack[-1])] = c_text
            self.emit( state, portion_barcode, portion_data, 'Portion' )
        """
        """
        for s_node, s_stack, s_attr, s_text in dom_scan(root_node, 'tcga_bcr/patient/samples/sample/portions/portion/analytes/analyte'):
            analyte_barcode = None
            for c_node, c_stack, c_attr, c_text in dom_scan(s_node, 'analyte/bcr_analyte_barcode'):
                analyte_barcode = c_text
            analyte_data = {}
            for c_node, c_stack, c_attr, c_text in dom_scan(s_node, 'analyte/*'):
                if 'xsd_ver' in c_attr:
                    analyte_data[c_attr.get('preferred_name', c_stack[-1])] = c_text
            self.emit( state, analyte_barcode, analyte_data, 'Analyte' )
        """
        """
        for s_node, s_stack, s_attr, s_text in dom_scan(root_node, 'tcga_bcr/patient/samples/sample/portions/portion/analytes/analyte/aliquots/aliquot'):
            aliquot_barcode = None
            for c_node, c_stack, c_attr, c_text in dom_scan(s_node, 'aliquot/bcr_aliquot_barcode'):
                aliquot_barcode = c_text
            aliquot_data = {}
            for c_node, c_stack, c_attr, c_text in dom_scan(s_node, 'aliquot/*'):
                if 'xsd_ver' in c_attr:
                    aliquot_data[c_attr.get('preferred_name', c_stack[-1])] = c_text
            self.emit( state, aliquot_barcode, aliquot_data, 'Aliquot' )
        """

        for s_node, s_stack, s_attr, s_text in dom_scan(root_node, 'tcga_bcr/patient/drugs/drug'):
            drug_barcode = None
            for c_node, c_stack, c_attr, c_text in dom_scan(s_node, 'drug/bcr_drug_barcode'):
                drug_barcode = c_text
            drug_data = deepcopy(patient_data)
            drug_data['source'] = source
            for c_node, c_stack, c_attr, c_text in dom_scan(s_node, 'drug/*'):
                if 'xsd_ver' in c_attr:
                    field_name = c_attr.get('preferred_name', c_stack[-1])
                    if len(field_name):
                        drug_data[field_name] = c_text
                    else:
                        drug_data[c_node.localName] = c_text
            self.emit( state, drug_barcode, drug_data, 'DrugTherapy' )

        for s_node, s_stack, s_attr, s_text in dom_scan(root_node, 'tcga_bcr/patient/radiations/radiation'):
            radiation_barcode = None
            for c_node, c_stack, c_attr, c_text in dom_scan(s_node, 'radiation/bcr_radiation_barcode'):
                radiation_barcode = c_text
            radiation_data = deepcopy(patient_data)
            radiation_data['source'] = source
            for c_node, c_stack, c_attr, c_text in dom_scan(s_node, 'radiation/*'):
                if 'xsd_ver' in c_attr:
                    field_name = c_attr.get('preferred_name', c_stack[-1])
                    if len(field_name):
                        radiation_data[field_name] = c_text
                    else:
                        radiation_data[c_node.localName] = c_text
            self.emit( state, radiation_barcode, radiation_data, 'RadiationTherapy' )

        for s_node, s_stack, s_attr, s_text in dom_scan(root_node, 'tcga_bcr/patient/follow_ups/follow_up'):
            follow_up_barcode = None
            sequence = s_attr['sequence']
            for c_node, c_stack, c_attr, c_text in dom_scan(s_node, 'follow_up/bcr_followup_barcode'):
                follow_up_barcode = c_text
            follow_up_data = deepcopy(patient_data)
            follow_up_data['sequence'] = sequence
            follow_up_data['source'] = source
            for c_node, c_stack, c_attr, c_text in dom_scan(s_node, 'follow_up/*'):
                if 'xsd_ver' in c_attr:
                    field_name = c_attr.get('preferred_name', c_stack[-1])
                    if len(field_name):
                        follow_up_data[field_name] = c_text
                    else:
                        follow_up_data[c_node.localName] = c_text
            self.emit( state, follow_up_barcode, follow_up_data, 'ClinicalFollowup' )

        return state

    def emit(self, state, key, entry, entry_type):
        state['generators'][entry_type].find(entry)

def extract_cohorts(state):
    for key in state['Individual']:
        individual = state['Individual'][key]
        project_code = individual.attributes['project_code']
        disease_code = individual.attributes['disease_code']
        cohort = state['generators']['Project'].find({'project_code': project_code, 'disease_code': disease_code})
        cohort.hasMember.append(individual.id)

def initial_state(pubchem={}):
    individual_generator = IndividualGenerator()
    biosample_generator = BiosampleGenerator(individual_generator.gid)
    cohort_generator = ProjectGenerator(individual_generator.gid)
    drug_generator = DrugGenerator(individual_generator.gid)
    drug_generator.pubchem = pubchem
    radiation_generator = RadiationGenerator(individual_generator.gid)
    followup_generator = FollowupGenerator(individual_generator.gid)

    return record_initial_state({
        'Individual': individual_generator,
        'Biosample': biosample_generator,
        'Project': cohort_generator,
        'DrugTherapy' : drug_generator,
        'RadiationTherapy' : radiation_generator,
        'ClinicalFollowup' : followup_generator
    })

def process_file(state, source, file):
    raw = file.read()
    try:
        print('parsing', file.name)
        dom = parseString(raw)
    except:
        dom = None

    if dom:
        extract = ClinicalParser()
        return extract.parseXMLFile(state, source, dom)
    else:
        return state


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--tar', type=str, default=None, help='a tarball clinical files')
    parser.add_argument('--input', type=str, default=None, help='path to directory containing clinical files')
    parser.add_argument('--output', type=str, help='path to output file')
    parser.add_argument('--source', type=str, help='source of data (eg. gdc)')
    parser.add_argument('--pubchem', type=str, help='pubchem id map')

    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = parse_arguments()

    pubchem_map = {}
    if args.pubchem:
        with open(args.pubchem) as handle:
            for line in handle:
                row = line.rstrip().split("\t")
                pubchem_map[row[0].lower()] = row[1]

    extract = ClinicalParser()
    state = initial_state(pubchem=pubchem_map)

    if args.input:
        process_input(state, args.input, args.source, process_file)
    elif args.tar:
        tmpdir = tempfile.mkdtemp(dir=".", prefix="gdc_extract_")
        subprocess.check_call("tar xvzf %s" % (os.path.abspath(args.tar)), shell=True, cwd=tmpdir)
        process_input(state, tmpdir, args.source, process_file)

    output_state(state, args.output)
