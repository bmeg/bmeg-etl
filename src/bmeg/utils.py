import inspect
import typing

from contextlib import suppress
from functools import wraps


def enforce_types(callable):
    spec = inspect.getfullargspec(callable)

    def check_types(*args, **kwargs):
        parameters = dict(zip(spec.args, args))
        parameters.update(kwargs)
        for name, value in parameters.items():
            # Assume un-annotated parameters can be any type
            with suppress(KeyError):
                type_hint = spec.annotations[name]
                if isinstance(type_hint, typing._SpecialForm):
                    # No check for typing.Any, typing.Union, typing.ClassVar
                    # (without parameters)
                    continue
                try:
                    actual_type = type_hint.__origin__
                except AttributeError:
                    actual_type = type_hint
                if isinstance(actual_type, typing._SpecialForm):
                    # case of typing.Union[…] or typing.ClassVar[…]
                    actual_type = type_hint.__args__

                if not isinstance(value, actual_type):
                    raise TypeError(
                        "Unexpected type for '{}' (expected {} but found {})".
                        format(name, type_hint, type(value))
                    )

    def decorate(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            check_types(*args, **kwargs)
            return func(*args, **kwargs)
        return wrapper

    if inspect.isclass(callable):
        callable.__init__ = decorate(callable.__init__)
        return callable

    return decorate(callable)


def set_gid(obj, gid):
    object.__setattr__(obj, "gid", gid)


def get_tcga_individual_barcode(id):
    parts = id.split("-")
    return "-".join(parts[0:3])


def get_tcga_sample_barcode(id):
    parts = id.split("-")
    return "-".join(parts[0:4])


def get_tcga_portion_barcode(id):
    parts = id.split("-")
    parts[5] = parts[5][:-1]
    return "-".join(parts[0:5])


def get_tcga_analyte_barcode(id):
    parts = id.split("-")
    return "-".join(parts[0:5])


def get_tcga_aliquot_barcode(id):
    parts = id.split("-")
    return "-".join(parts[0:7])


def tcga_barcode_is_tumor(id):
    parts = id.split("-")
    sample_number = parts[4][:-1]
    return sample_number < 10


def tcga_barcode_is_normal(id):
    parts = id.split("-")
    sample_number = parts[4][:-1]
    return sample_number >= 10
