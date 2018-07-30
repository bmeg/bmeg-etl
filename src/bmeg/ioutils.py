import csv
import gzip
import io


def reader(path):
    if path.endswith(".gz"):
        return io.TextIOWrapper(gzip.GzipFile(path))
    else:
        return open(path, "r")


def csv(path, **kwargs):
    r = reader(path)
    return csv.DictReader(r, **kwargs)


def tsv(path, **kwargs):
    r = reader(path)
    return csv.DictReader(r, delimiter="\t", **kwargs)
