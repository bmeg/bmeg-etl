def assign_attributes(obj, init_vars):
    for name in init_vars:
        object.__setattr__(obj, name, init_vars[name])


def set_gid(obj, gid):
    object.__setattr__(obj, "gid", gid)
