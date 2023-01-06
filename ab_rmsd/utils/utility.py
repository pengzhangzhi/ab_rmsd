def exists(x):
    return x is not None


def exist_key(d, k):
    return d.get(k) is not None

class PDBParseError(Exception):
    pass