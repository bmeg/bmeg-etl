import dataclasses
import hashlib
import jsonschema
import logging
import os
import pkg_resources
import re
import types

from copy import deepcopy
from dictionaryutils import DataDictionary, load_schemas_from_dir, load_yaml
from functools import partial

from bmeg.gid import gid_factories, default_gid, cast_gid
from bmeg.utils import enforce_types


class BMEGDataDictionary(DataDictionary):
    """
    Modified from:
    https://github.com/uc-cdis/dictionaryutils/blob/42bf330d82bf084141c0f21b9815cc7e34bf5287/dictionaryutils/__init__.py#L112
    """

    def __init__(
            self,
            root_dir,
            definitions_paths=None,
            metaschema_path=None
    ):
        self.root_dir = root_dir
        self.metaschema_path = metaschema_path or self._metaschema_path
        self.definitions_paths = definitions_paths or self._definitions_paths
        self.exclude = (
            [self.metaschema_path] + self.definitions_paths + [self.settings_path] + ["data_release.yaml", "root.yaml"]
        )
        self.schema = dict()
        self.resolvers = dict()

        self.metaschema = load_yaml(
            os.path.join(
                pkg_resources.resource_filename("dictionaryutils", "schemas"),
                self.metaschema_path
            )
        )

        self.load_data(directory=self.root_dir, url=None)

    def load_data(self, directory=None, url=None):
        """Load and reslove all schemas from directory or url"""
        yamls, resolvers = load_schemas_from_dir(pkg_resources.resource_filename("dictionaryutils", "schemas/"))
        yamls, resolvers = load_schemas_from_dir(directory,
                                                 schemas=yamls,
                                                 resolvers=resolvers)

        self.settings = yamls.get(self.settings_path) or {}
        self.resolvers.update(resolvers)

        schemas = {
            schema["id"]: self.resolve_schema(schema, deepcopy(schema))
            for path, schema in yamls.items()
            if path not in self.exclude
        }
        self.schema.update(schemas)


def validate(props, schema):
    jsonschema.validate(props, schema)
    return


def get_props(in_data, preserve_null=False, keep=[], remove=[]):
    data = in_data.copy()
    to_keep = keep.copy()
    to_remove = remove.copy()
    if not preserve_null:
        nulls = [k for k in data if data[k] is None and k not in to_remove and k not in to_keep]
        to_remove += nulls
        for k in to_remove:
            del data[k]
    return data


class ClassInstance:

    def __init__(self, **kwargs):
        self.__dict__["_props"] = {}
        if "properties" in self._schema:
            for k in self._schema["properties"].keys():
                self.__dict__["_props"][k] = None

        for k, v in kwargs.items():
            self.__setattr__(k, v)

    def props(self, preserve_null=False, keep=["submitter_id"], remove=["id"]):
        return get_props(in_data=self._props, preserve_null=preserve_null, keep=keep, remove=remove)

    def schema(self):
        return self._schema

    def validate(self):
        validate(self.props(), self.schema())
        return

    def label(self):
        return self.__class__.__name__

    def gid(self):
        if not isinstance(self.id, self._gid_cls):
            raise TypeError("expected id of type {}", self._gid_cls)
        return self.id

    def __repr__(self):
        return '<%s(%s)>' % (self.__class__.__name__,
                             self.props(preserve_null=True))

    def __setattr__(self, key, item):
        if key not in self._props:
            raise KeyError("object does not contain key '{}'".format(key))
        self._props[key] = item

    def __getattr__(self, key):
        return self._props[key]

    def __getitem__(self, k):
        return self.__getattr__(k)

    def __setitem__(self, k, v):
        return self.__setattr__(k, v)


def capitalize(label):
    label = "".join([(lambda x: x[0].upper() + x[1:])(x) for x in re.split("_| +", label)])
    if label.startswith("G2p"):
        label = "G2P" + label[3:]
    return label


class Vertex:
    pass


class Edge:
    pass


@enforce_types
@dataclasses.dataclass(frozen=True)
class Deadletter(Vertex):
    """ standard way to log missing data """
    target_label: str  # desired vertex label
    data: dict  # data from source

    def validate(self):
        return

    def label(self):
        return self.__class__.__name__

    def gid(self):
        return Deadletter.make_gid(self.target_label, self.data)

    def props(self):
        return self.__dict__

    @classmethod
    def make_gid(cls, target_label, data):
        # create hash from data dict
        data = '%s:%s' % (target_label, str(data))
        datahash = hashlib.sha1()
        datahash.update(data.encode('utf-8'))
        datahash = datahash.hexdigest()
        return "%s:%s" % (cls.__name__, datahash)


_schemaPath = pkg_resources.resource_filename(__name__, "bmeg-dictionary/gdcdictionary/schemas")
_schema = BMEGDataDictionary(root_dir=_schemaPath)

__all__ = ['Vertex', 'Edge', 'Deadletter']
for k, schema in _schema.schema.items():
    name = capitalize(schema["id"])
    if name not in gid_factories:
        logging.debug("using default gid factory for '{}'".format(name))
        gid_factory = partial(default_gid, name)
    else:
        gid_factory = gid_factories[name]
    cls = type(
        name,
        (ClassInstance, Vertex),
        {
            '_schema': schema,
            '_gid_cls': types.new_class("{}GID".format(name), (str,), {}),
            'make_gid': classmethod(cast_gid(gid_factory))
        }
    )
    globals()[name] = cls
    __all__.append(name)

for k, schema in _schema.schema.items():
    for link in schema['links']:
        # TODO: handle link subgroup?
        if "target_type" not in link:
            continue
        # should we add required fields?
        edge_schema = {
            '$schema': 'http://json-schema.org/draft-04/schema#',
            'type': 'object',
            'additionalProperties': False,
            'properties': link.get('properties', {})
        }
        src = capitalize(k)
        target = capitalize(link["target_type"])
        cls_name = "{}_{}_{}".format(src, capitalize(link['label']), target)
        cls = enforce_types(
            dataclasses.make_dataclass(
                cls_name=cls_name,
                fields=[
                    ('from_gid', globals()[src]._gid_cls),
                    ('to_gid', globals()[target]._gid_cls),
                    ('__label', str, dataclasses.field(default=link['label'], init=False)),
                    ('data', dict, dataclasses.field(default_factory=dict))
                ],
                bases=(Edge,),
                namespace={
                    '_schema': edge_schema,
                    'schema': lambda self: self._schema,
                    'validate': lambda self: validate(self.props(), self.schema()),
                    'label': lambda self: self.__label,
                    'props': lambda self: get_props(self.data)
                }
            )
        )
        if 'backref' not in link:
            cls._backref = None
            cls.backref = lambda self: None
            globals()[cls_name] = cls
            __all__.append(cls_name)
            continue
        bkref = "{}_{}_{}".format(target, capitalize(link['backref']), src)
        bkref_cls = enforce_types(
            dataclasses.make_dataclass(
                cls_name=bkref,
                fields=[
                    ('from_gid', globals()[target]._gid_cls),
                    ('to_gid', globals()[src]._gid_cls),
                    ('__label', str, dataclasses.field(default=link['backref'], init=False)),
                    ('data', dict, dataclasses.field(default_factory=dict))
                ],
                bases=(Edge,),
                namespace={
                    '_schema': edge_schema,
                    'schema': lambda self: self._schema,
                    'validate': lambda self: validate(self.props(), self.schema()),
                    'label': lambda self: self.__label,
                    'props': lambda self: get_props(self.data)
                }
            )
        )
        cls._backref = bkref_cls
        cls.backref = lambda self: self._backref(from_gid=self.to_gid, to_gid=self.from_gid)
        bkref_cls._backref = cls
        bkref_cls.backref = lambda self: self._backref(from_gid=self.to_gid, to_gid=self.from_gid)
        globals()[cls_name] = cls
        globals()[bkref] = bkref_cls
        __all__.append(cls_name)
        __all__.append(bkref)
