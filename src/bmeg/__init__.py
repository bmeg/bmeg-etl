import pkg_resources
import jsonschema
import os
import types

from dictionaryutils import DataDictionary, load_schemas_from_dir, load_yaml

from copy import deepcopy


class BMEGDataDictionary(DataDictionary):
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
            [self.metaschema_path] +
            self.definitions_paths +
            [self.settings_path] +
            ["data_release.yaml", "root.yaml"]
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


class ClassInstance:

    def __init__(self, **kwargs):
        self.__dict__["_props"] = {}
        if "properties" in self._schema:
            for k in self._schema["properties"].keys():
                self.__dict__["_props"][k] = None

        for k, v in kwargs.items():
            self.__setattr__(k, v)

    def props(self):
        return self._props

    def schema(self):
        return self._schema

    def validate(self):
        jsonschema.validate(self.props(), self.schema())
        return

    def label(self):
        return self.__class__.__name__

    def gid(self):
        return ClassInstance.make_gid(self.id)

    @classmethod
    def make_gid(cls, gid):
        return cls._gid_cls(gid)

    def __repr__(self):
        return '<%s(%s)>' % (self.__class__.__name__, self.props())

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


_schemaPath = pkg_resources.resource_filename(__name__, "bmeg-dictionary/gdcdictionary/schemas")
_schema = BMEGDataDictionary(root_dir=_schemaPath)

__all__ = []
for k, schema in _schema.schema.items():
    name = k.capitalize()
    cls = type(
        name, (ClassInstance,),
        {'_schema': schema,
         '_gid_cls': types.new_class("{}GID".format(name), (str,), {})}
    )
    globals()[name] = cls
    __all__.append(name)
