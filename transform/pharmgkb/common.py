import stringcase


def snake_case_keys(line):
    """Transforms sentence case to snake case."""
    keys = []
    for k in line.keys():
        _k = stringcase.snakecase(k).replace('__', '_')
        keys.append(_k)
    return keys
