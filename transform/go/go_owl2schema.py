#!/usr/bin/env python

# http://purl.obolibrary.org/obo/go/extensions/go-plus.owl

import sys
import re
import subprocess
from bmeg import phenotype_pb2

def which(cmd):
    cmd = ["which",cmd]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    res = p.stdout.readline().rstrip()
    if len(res) == 0: return None
    return res


re_string_type = re.compile(r'"(.*)"\^\^<http://www.w3.org/2001/XMLSchema\#string>')
re_bool_type = re.compile(r'"(.*)"\^\^<http://www.w3.org/2001/XMLSchema\#boolean>')

def type_convert(s):
    res = re_string_type.search(s)
    if res:
        return res.group(1)
    res = re_bool_type.search(s)
    if res:
        return bool(res.group(1))
    return None

def follow_intersection(elem, objects):
    if 'intersectionOf' in elem:
        o = {'class' : 'Intersection'}
        n = objects[ elem['intersectionOf'][0] ]
        members, visited = follow_intersection(n, objects)
        o['memberEdges'] = members
        return o, visited
    elif 'first' in elem:
        
        o = [ elem['first'][0] ]
        
        visited = []
        if 'rest' in elem and elem['rest'][0] in objects:
             n = objects[elem['rest'][0]]
             members, visited = follow_intersection(n, objects)
             o += members
        return o, visited + [ elem['rest'][0] ]
    else:
        print elem
        raise Exception("unknown")

if __name__ == "__main__":
    objects = {}
    proc = subprocess.Popen([which("rapper"), sys.argv[1]], stdout=subprocess.PIPE)
    re_hash_link = re.compile(r'<.*\#(.*)>')
    re_slash_link = re.compile(r'<.*\/(.*)>')
    re_turtle_end = re.compile(r' \.\s*$')
    for i in proc.stdout:
        #break every line into subject, predicate, object
        tmp = re_turtle_end.sub("", i).split(" ")
        sub = tmp[0]
        pred = tmp[1]
        obj = " ".join(tmp[2:])
        
        res = re_hash_link.search(pred)
        if res is None:
            res = re_slash_link.search(pred)
        
        pred_id = res.group(1)
        
        if sub not in objects:
            objects[sub] = {}
        if pred_id not in objects[sub]:
            objects[sub][pred_id] = []
        
        v = type_convert(obj)
        if v is not None:
            objects[sub][pred_id].append(v)
        else:
            objects[sub][pred_id].append(obj)
    proc.communicate()
    
    removed = set()
    outset = {}
    for rec_id, rec in objects.items():
        if rec_id not in removed:
            o = phenotype_pb2.OntologyTerm()
            if 'id' in rec:
                o.gid = rec['id'][0]
            if 'type' in rec:
                oClass = re_hash_link.search(rec['type'][0]).group(1)
                if oClass == 'Class' and 'intersectionOf' in rec:
                    isect, visited = follow_intersection(rec, objects)
                    outset[rec_id] = isect
                    for i in visited:
                        removed.add(i)
                    o = None
            if o is not None:
                for k, v in rec.items():
                    if k not in ['id', 'type']:
                        o[k] = v
                outset[rec_id] = o
    for rec_id, rec in outset.items():
        if rec_id not in removed:
            print rec_id, rec