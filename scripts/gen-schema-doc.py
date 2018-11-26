#!/usr/bin/env python

import re
import inspect

from bmeg import vertex
from bmeg import edge

print("# Vertices")
for v in dir(vertex):
    c = getattr(vertex, v)
    if isinstance(c, type):
        if issubclass(c, vertex.Vertex):
            name = v
            #print(dir(c))
            desc = inspect.getcomments(c)
            if desc is None:
                desc = "N/A"
            else:
                desc = re.sub(r'^\#', '', desc).rstrip()
            print("## %s : %s" % (name, desc))
            for k, val in c.__dataclass_fields__.items():
                typename = "None"
                try:
                    typename = val.type.__name__
                except:
                    pass
                desc = inspect.getcomments(val) #BUG: This doesn't seem to work on Fields
                if desc is None:
                    desc = "N/A"
                print(" - %s (%s) : %s " % (k, typename, desc))
            print()

print("# Edges")
for e in dir(edge):
    c = getattr(edge, e)
    if isinstance(c, type):
        if issubclass(c, edge.Edge):
            name = e
            #print(dir(c))
            print("##" + name)
            for k, v in c.__dataclass_fields__.items():
                #print("  " + k)
                #try:
                #print(v)
                #pdb.set_trace()
                #print(dir(v.type))
                typename = "None"
                try:
                    typename = e.type.__name__
                except:
                    pass
                print(" - %s (%s)" % (k, typename))
            print()
