import csv
import gzip
import io


def reader(path):
    if path.endswith(".gz"):
        return io.TextIOWrapper(gzip.GzipFile(path))
    else:
        return open(path, "r")


def read_csv(path, **kwargs):
    r = reader(path)
    return csv.DictReader(r, **kwargs)


def read_tsv(path, **kwargs):
    r = reader(path)
    return csv.DictReader(r, delimiter="\t", **kwargs)
