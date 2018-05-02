import os
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

def process_clinical(state, input):
    raw = input.read()

    try:
        dom = parseString(raw)
        root = dom.childNodes[0]
        code = 'TCGA'
        site = 'body'

        for node, stack, attr, text in dom_scan(root, 'tcga_bcr/admin/disease_code'):
            code = text

        for node, stack, attr, text in dom_scan(root, 'tcga_bcr/patient/tumor_tissue_site'):
            site = text

        if code in state['mapping']:
            if not site == state['mapping'][code]:
                print('differing sites for project code!', code, site)
        else:
            state['mapping'][code] = site

        if code in state['diffs']:
            state['diffs'][code].add(site)
        else:
            state['diffs'][code] = set([site])

        print(code, site)
    except:
        print('failed to parse', input.name)

    return state

def process_input(path, process):
    state = {'mapping': {}, 'diffs': {}}
    if os.path.isdir(path):
        for item in os.listdir(path):
            with open(os.path.join(path, item)) as input:
                state = process(state, input)

    for key in state['diffs'].keys():
        state['diffs'][key] = list(state['diffs'][key])

    print json.dumps(state)

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=str, help='path to file containing clinical attributes')
    parser.add_argument('--output', type=str, help='path to output file')
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = parse_arguments()
    process_input(args.input, process_clinical)
