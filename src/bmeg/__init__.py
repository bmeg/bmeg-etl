
import os
import yaml
import json
from copy import copy
from glob import glob



def load(path):
    out = Schema()
    for g in glob(os.path.join(path, "*.yaml")):
        with open(g) as handle:
            meta = yaml.load(handle.read(), Loader=yaml.SafeLoader)
        if meta.get("type", "") == "object":
            c = SchemaClass(out, os.path.basename(g), meta)
        else:
            out._add(os.path.basename(g), meta)
    return out

class Schema(object):
    def __init__(self):
        self._classes = {}
        self._path = {}

    def _add(self, path, c):
        if isinstance(c, SchemaClass):
            self._classes[c.name] = c
        self._path[path] = c

    def _getTerm(self, p):
        f, elem = p.split("#")
        if f in self._path:
            meta = self._path[f]
            elem = elem.lstrip("/")
            if elem in meta:
                return meta[elem]
        return None

    def _getProp(self, p):
        f, elem = p.split("#")
        if f in self._path:
            meta = self._path[f]
            elem = elem.lstrip("/")
            if elem in meta:
                return meta[elem]

    def __getattr__(self,name):
        if name in self._classes:
            return self._classes[name]
        raise AttributeError("%s" % (name))

class SchemaClass:

    def __init__(self, schema, path, classDict):
        self._schema = schema
        self._path = path
        self._classDict = classDict

        self._schema._add(path, self)

    def __call__(self, **data):
        return ClassInstance(self._schema, self, **data)

    def getName(self):
        #print(self._schema)
        return self._classDict['title']

    def getLink(self, n):
        for i in self._classDict.get("links", []):
            if i['name'] == n:
                return i
        return None

    def properties(self):
        return self._classDict['properties']

    def prop(self, name):
        prop = self._classDict['properties'][name]
        if "$ref" in prop:
            m = copy(self._schema._getProp(prop["$ref"]))
            for k, v in prop.items():
                if k != "$ref":
                    m[k] = v
            return SchemaProperty(self._schema, m)
        else:
            return SchemaProperty(self._schema, prop)

    def required(self):
        return self._classDict["required"]

    name = property(getName)
    label = property(getName)

class ClassInstance:
    def __init__(self, schema, classSchema, **data):
        self._schema = schema
        self._classSchema = classSchema
        self._data = {}

        for k, v in data.items():
            o = self._validateSet(k, v)
            if o is not None:
                print(o)

    def _validateSet(self, k, v):
        if k not in self._classSchema.properties():
            return "Property '%s' not found" % (k)
        prop = self._classSchema.prop(k)
        vt = None
        ve = []
        for t in prop.types:
            enum = t.enum
            if len(enum) > 0:
                if v not in enum:
                    ve.append( "value %s not in [%s]" % (v, ",".join(enum)) )
                else:
                    vt = t
            else:
                vt = t
        if vt is not None:
            self._data[k] = v
        else:
            return "type not found: [%s]" % (",".join(ve))

    def _validate(self):
        o = []
        for i in self._classSchema.required():
            if i not in self._data:
                o.append("required property %s not found" % (i))
        return o

    def __getattr__(self,name):
        if name in self._data:
            return self._data[name]
        raise AttributeError("%s" % (name))

class SchemaProperty:
    def __init__(self, schema, propDict):
        self._schema = schema
        self._propDict = propDict

    def getSystemAlias(self):
        #print(self._propDict)
        return self._propDict.get("systemAlias", "")

    #def types(self):
    #return self._propDict.get("oneOf", [])

    def getTypes(self):
        if "oneOf" in self._propDict:
            o = []
            for i in self._propDict["oneOf"]:
                o.append(SchemaType(self._schema, i))
            return o
        return [ SchemaType(self._schema, self._propDict) ]

    systemAlias = property(getSystemAlias)
    types = property(getTypes)

class SchemaType:
    def __init__(self, schema, typeDict):
        self._schema = schema
        self._typeDict = typeDict

    def getEnum(self):
        return self._typeDict.get("enum", [])

    enum = property(getEnum)

_schema = load( os.path.join(os.path.dirname(os.path.abspath(__file__)), "schemas") )

__all__ = []
for k, v in _schema._classes.items():
    globals()[k] = v
    __all__.append(k)
