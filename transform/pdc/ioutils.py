"""Useful io utilities."""

import json
import gzip
import csv
import io


class JsonReader:
    def __init__(self, path):
        if path.endswith(".json.gz"):
            self.fp = io.TextIOWrapper(io.BufferedReader(gzip.GzipFile(path)))
        else:
            self.fp = open(path, "r", encoding='utf-8')

    def __iter__(self):
        return self

    def __next__(self):
        try:
            line = self.fp.readline()
            if len(line) < 1:
                raise IndexError()
            return json.loads(line)
        except IndexError:
            raise StopIteration()

    def __enter__(self):
        # create logic, already opened
        return self

    def __exit__(self, type, value, traceback):
        # ensure file closed
        self.fp.close()


class GzipBufferedTextReader(io.TextIOWrapper):
    """Wraps TextIOWrapper with buffered GzipFile reader."""
    def __init__(self, file, **kwargs):
        super().__init__(io.BufferedReader(gzip.GzipFile(file)), **kwargs)


class CsvReader(csv.DictReader):
    """Wraps csv.DictReader with tsv options."""
    def __init__(self, file, **kwargs):
        self.fp = open(file, "r")
        if file.endswith('.tsv'):
            super().__init__(self.fp, delimiter="\t", **kwargs)
        else:
            super().__init__(self.fp, **kwargs)

    def __enter__(self):
        # create logic, already opened
        return self

    def __exit__(self, type, value, traceback):
        # ensure file closed
        self.fp.close()


class GzipBufferedDictReader(csv.DictReader):
    """Wraps csv.DictReader with buffered GzipFile reader."""
    def __init__(self, file, **kwargs):
        if '.tsv' in file and 'delimiter' not in kwargs:
            kwargs['delimiter'] = '\t'
        self.fp = GzipBufferedTextReader(file)
        super().__init__(self.fp, **kwargs)

    def __enter__(self):
        # create logic, already opened
        return self

    def __exit__(self, type, value, traceback):
        # ensure file closed
        self.fp.close()


def recommended_reader(path):
    """Deduces a recommended reader class based on path."""
    if path.endswith(".csv.gz") or path.endswith(".tsv.gz"):
        return GzipBufferedDictReader
    if path.endswith(".json.gz") or path.endswith(".json"):
        return JsonReader
    if path.endswith(".gz"):
        return GzipBufferedTextReader
    if path.endswith(".tsv") or path.endswith(".csv"):
        return CsvReader
    return None


def reader(path, reader_class=None, **kwargs):
    """Construct a file reader object from file specified by path,
    or a DictReader object if the file is a .tsv, .csv, .json.
    Automatically handles gz
    """
    # no preference for reader specified, look up one of ours
    if not reader_class:
        reader_class = recommended_reader(path)
    # no reader, open in text mode
    if not reader_class:
        return open(path, "r")
    return reader_class(path, **kwargs)
