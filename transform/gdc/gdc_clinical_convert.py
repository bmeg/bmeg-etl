#!/usr/bin/env python

import json
import argparse
from xml.dom.minidom import parseString


def getText(nodelist):
    rc = []
    for node in nodelist:
        if node.nodeType == node.TEXT_NODE:
            rc.append(node.data)
    return ''.join(rc)

def dom_scan(node, query):
    stack = query.split("/")
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

class ClinicalParser:

    def __init__(self):
        pass

    def parseXMLFile(self, dom, dataSubType):
        root_node = dom.childNodes[0]
        admin = {}
        for node, stack, attr, text in dom_scan(root_node, "tcga_bcr/admin/*"):
            admin[stack[-1]] = text

        patient_barcode = None
        for node, stack, attr, text in dom_scan(root_node, 'tcga_bcr/patient/bcr_patient_barcode'):
            patient_barcode = text

        patient_data = {}
        for node, stack, attr, text in dom_scan(root_node, "tcga_bcr/patient/*"):
            if 'xsd_ver' in attr:
                #print patientName, stack[-1], attr, text
                p_name = attr.get('preferred_name', stack[-1])
                if len(p_name) == 0:
                    p_name = stack[-1]
                patient_data[p_name] = text

        #if dataSubType == "patient":
        for node, stack, attr, text in dom_scan(root_node, "tcga_bcr/patient/stage_event/*"):
            if 'xsd_ver' in attr:
                p_name = attr.get('preferred_name', stack[-1])
                if len(p_name) == 0:
                    p_name = stack[-1]
                patient_data[p_name] = text
        for node, stack, attr, text in dom_scan(root_node, "tcga_bcr/patient/stage_event/*/*"):
            if 'xsd_ver' in attr:
                p_name = attr.get('preferred_name', stack[-1])
                if len(p_name) == 0:
                    p_name = stack[-1]
                patient_data[p_name] = text
        for node, stack, attr, text in dom_scan(root_node, "tcga_bcr/patient/stage_event/tnm_categories/*/*"):
            if 'xsd_ver' in attr:
                p_name = attr.get('preferred_name', stack[-1])
                if len(p_name) == 0:
                    p_name = stack[-1]
                patient_data[p_name] = text
        self.emit( patient_barcode, patient_data, "Patient" )

        #if dataSubType == "sample":
        for s_node, s_stack, s_attr, s_text in dom_scan(root_node, "tcga_bcr/patient/samples/sample"):
            sample_barcode = None
            for c_node, c_stack, c_attr, c_text in dom_scan(s_node, "sample/bcr_sample_barcode"):
                sample_barcode = c_text
            sample_data = {}
            for c_node, c_stack, c_attr, c_text in dom_scan(s_node, "sample/*"):
                if 'xsd_ver' in c_attr:
                    sample_data[c_attr.get('preferred_name', c_stack[-1])] = c_text
            self.emit( sample_barcode, sample_data, "Sample" )

        #if dataSubType == "portion":
        for s_node, s_stack, s_attr, s_text in dom_scan(root_node, "tcga_bcr/patient/samples/sample/portions/portion"):
            portion_barcode = None
            for c_node, c_stack, c_attr, c_text in dom_scan(s_node, "portion/bcr_portion_barcode"):
                portion_barcode = c_text
            portion_data = {}
            for c_node, c_stack, c_attr, c_text in dom_scan(s_node, "portion/*"):
                if 'xsd_ver' in c_attr:
                    portion_data[c_attr.get('preferred_name', c_stack[-1])] = c_text
            self.emit( portion_barcode, portion_data, "Portion" )

        #if dataSubType == "analyte":
        for s_node, s_stack, s_attr, s_text in dom_scan(root_node, "tcga_bcr/patient/samples/sample/portions/portion/analytes/analyte"):
            analyte_barcode = None
            for c_node, c_stack, c_attr, c_text in dom_scan(s_node, "analyte/bcr_analyte_barcode"):
                analyte_barcode = c_text
            analyte_data = {}
            for c_node, c_stack, c_attr, c_text in dom_scan(s_node, "analyte/*"):
                if 'xsd_ver' in c_attr:
                    analyte_data[c_attr.get('preferred_name', c_stack[-1])] = c_text
            self.emit( analyte_barcode, analyte_data, "Analyte" )

        #if dataSubType == "aliquot":
        for s_node, s_stack, s_attr, s_text in dom_scan(root_node, "tcga_bcr/patient/samples/sample/portions/portion/analytes/analyte/aliquots/aliquot"):
            aliquot_barcode = None
            for c_node, c_stack, c_attr, c_text in dom_scan(s_node, "aliquot/bcr_aliquot_barcode"):
                aliquot_barcode = c_text
            aliquot_data = {}
            for c_node, c_stack, c_attr, c_text in dom_scan(s_node, "aliquot/*"):
                if 'xsd_ver' in c_attr:
                    aliquot_data[c_attr.get('preferred_name', c_stack[-1])] = c_text
            self.emit( aliquot_barcode, aliquot_data, "Aliquot" )

        #if dataSubType == "drug":
        for s_node, s_stack, s_attr, s_text in dom_scan(root_node, "tcga_bcr/patient/drugs/drug"):
            drug_barcode = None
            for c_node, c_stack, c_attr, c_text in dom_scan(s_node, "drug/bcr_drug_barcode"):
                drug_barcode = c_text
            drug_data = {}
            for c_node, c_stack, c_attr, c_text in dom_scan(s_node, "drug/*"):
                if 'xsd_ver' in c_attr:
                    drug_data[c_attr.get('preferred_name', c_stack[-1])] = c_text
            self.emit( drug_barcode, drug_data, "Drug" )

        #if dataSubType == "radiation":
        for s_node, s_stack, s_attr, s_text in dom_scan(root_node, "tcga_bcr/patient/radiations/radiation"):
            radiation_barcode = None
            for c_node, c_stack, c_attr, c_text in dom_scan(s_node, "radiation/bcr_radiation_barcode"):
                radiation_barcode = c_text
            radiation_data = {}
            for c_node, c_stack, c_attr, c_text in dom_scan(s_node, "radiation/*"):
                if 'xsd_ver' in c_attr:
                    radiation_data[c_attr.get('preferred_name', c_stack[-1])] = c_text
            self.emit( radiation_barcode, radiation_data, "Radiation" )

        #if dataSubType == "followup":
        for s_node, s_stack, s_attr, s_text in dom_scan(root_node, "tcga_bcr/patient/follow_ups/follow_up"):
            follow_up_barcode = None
            sequence = s_attr['sequence']
            for c_node, c_stack, c_attr, c_text in dom_scan(s_node, "follow_up/bcr_followup_barcode"):
                follow_up_barcode = c_text
            follow_up_data = { "sequence" : sequence}
            for c_node, c_stack, c_attr, c_text in dom_scan(s_node, "follow_up/*"):
                if 'xsd_ver' in c_attr:
                    follow_up_data[c_attr.get('preferred_name', c_stack[-1])] = c_text
            self.emit( follow_up_barcode, follow_up_data, "Followup" )
    def emit(self, key, entry, entryType):
        out = {
            "gid" : "tcga:%s" % key,
            "type" : entryType,
            "meta" : {}
        }
        for k, v in entry.items():
            if len(v):
                out['meta'][k] = v
        print json.dumps(out)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("subtype")
    parser.add_argument("file")

    args = parser.parse_args()

    clin = ClinicalParser()
    with open(args.file) as handle:
        data = handle.read()
    dom = parseString(data)
    clin.parseXMLFile(dom, args.subtype)
