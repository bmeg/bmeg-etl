#!/usr/bin/env python3

import argparse
import gzip
import json
import logging
import re
import sys
from glob import glob
import xml.sax
import os
from ftplib import FTP
from multiprocessing import Pool, Queue, Manager
from functools import partial

from bmeg import Publication, Project
from bmeg.emitter import JSONEmitter


reWord = re.compile(r'\w')
reSpace = re.compile(r'\s')


def ignore(e, v, attrs, **kwds):
    return None


def create_list(e, v, attrs):
    return [v]


def string_pass(e, v, attrs):
    return v


def pass_data(e, v, attrs, **kwds):
    return kwds


def pass_attrs(e, v, attrs, **kwds):
    attrs.update(kwds)
    return attrs


def create_dict_list(e, v, attrs, **kwds):
    return [kwds]


def date_extract(e, v, attrs, **kwds):
    # if 'MedlineDate' in kwds:
    #     return kwds['MedlineDate']
    if 'Year' in kwds and 'Month' in kwds and 'Day' in kwds:
        return "%s-%s-%s" % (kwds['Year'], kwds['Month'], kwds['Day'])
    elif 'Year' in kwds and 'Month' in kwds:
        return "%s-%s" % (kwds['Year'], kwds['Month'])
    elif 'Year' in kwds:
        return kwds['Year']
    else:
        return ''


def author_extract(e, v, attrs, **kwds):
    if 'ForeName' in kwds and 'LastName' in kwds:
        return ["%s %s" % (kwds['ForeName'], kwds['LastName'])]
    elif 'CollectiveName' in kwds:
        return [kwds['CollectiveName']]


def emit_pubmed(e, v, attrs, **kwds):
    # https://www.nlm.nih.gov/bsd/licensee/elements_descriptions.html#status_value
    pmid = kwds['MedlineCitation']['PMID']
    title = kwds['MedlineCitation']['Article']['ArticleTitle']
    date = kwds['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate']

    author = []
    if 'AuthorList' in kwds['MedlineCitation']['Article']:
        if 'Author' in kwds['MedlineCitation']['Article']['AuthorList']:
            [author.append(x) for x in kwds['MedlineCitation']['Article']['AuthorList']['Author']]
            if kwds['MedlineCitation']['Article']['AuthorList']['CompleteYN'] == 'N':
                author.append('et al.')
    abstract = ""
    if 'Abstract' in kwds['MedlineCitation']['Article']:
        abstract = kwds['MedlineCitation']['Article']['Abstract']['AbstractText']
    url = 'https://www.ncbi.nlm.nih.gov/pubmed/{}'.format(pmid)
    #out = Publication(id=Publication.make_gid(url),
    #                  url=url, title=title, abstract=abstract,
    #                  text="", date=date, author=author, citation=[],
    #                  project_id=Project.make_gid("Reference"))
    out = {
        "url":url,
        "title":title, "abstract":abstract,
        "text":"", "date":date, "author":author, "citation":[]
    }
    e.append(out)


f_map = [
    (['PubmedArticleSet', 'PubmedArticle'], None, emit_pubmed),
    (['PubmedArticleSet', 'PubmedArticle', 'MedlineCitation'], None, pass_attrs),
    (['PubmedArticleSet', 'PubmedArticle', 'MedlineCitation', 'PMID'], None, string_pass),
    (['PubmedArticleSet', 'PubmedArticle', 'MedlineCitation', 'DateCreated', '*'], None, string_pass),
    (['PubmedArticleSet', 'PubmedArticle', 'MedlineCitation', 'DateCreated'], None, date_extract),
    (['PubmedArticleSet', 'PubmedArticle', 'MedlineCitation', 'DateCompleted', '*'], None, string_pass),
    (['PubmedArticleSet', 'PubmedArticle', 'MedlineCitation', 'DateCompleted'], None, date_extract),
    (['PubmedArticleSet', 'PubmedArticle', 'MedlineCitation', 'DateRevised', '*'], None, string_pass),
    (['PubmedArticleSet', 'PubmedArticle', 'MedlineCitation', 'DateRevised'], None, date_extract),
    (['PubmedArticleSet', 'PubmedArticle', 'MedlineCitation', 'Article'], None, pass_data),
    (['PubmedArticleSet', 'PubmedArticle', 'MedlineCitation', 'Article', 'ArticleTitle'], None, string_pass),
    (['PubmedArticleSet', 'PubmedArticle', 'MedlineCitation', 'Article', 'Abstract'], None, pass_data),
    (['PubmedArticleSet', 'PubmedArticle', 'MedlineCitation', 'Article', 'Abstract', 'AbstractText'], None, string_pass),
    (['PubmedArticleSet', 'PubmedArticle', 'MedlineCitation', 'Article', 'AuthorList'], None, pass_attrs),
    (['PubmedArticleSet', 'PubmedArticle', 'MedlineCitation', 'Article', 'AuthorList', 'Author', '*'], None, string_pass),
    (['PubmedArticleSet', 'PubmedArticle', 'MedlineCitation', 'Article', 'AuthorList', 'Author'], None, author_extract),
    (['PubmedArticleSet', 'PubmedArticle', 'MedlineCitation', 'Article', 'Journal'], None, pass_data),
    (['PubmedArticleSet', 'PubmedArticle', 'MedlineCitation', 'Article', 'Journal', 'JournalIssue'], None, pass_data),
    (['PubmedArticleSet', 'PubmedArticle', 'MedlineCitation', 'Article', 'Journal', 'JournalIssue', 'PubDate', '*'], None, string_pass),
    (['PubmedArticleSet', 'PubmedArticle', 'MedlineCitation', 'Article', 'Journal', 'JournalIssue', 'PubDate'], None, date_extract),
]


class StackLevel:
    def __init__(self, name, attrs):
        self.data = {}
        self.name = name
        self.has_children = False
        self.attrs = attrs


def stack_match(query, elem):
    if len(query) != len(elem):
        return False
    for q, e in zip(query, elem):
        if isinstance(q, list):
            if e not in q:
                return False
        else:
            if q != "*" and q != e:
                return False
    return True


NOT_FOUND = set()


class PubMedHandler(xml.sax.ContentHandler):
    def __init__(self, record_write):
        xml.sax.ContentHandler.__init__(self)
        self.record_write = record_write
        self.stack = []

    def startElement(self, name, attrs):
        if len(self.stack):
            self.stack[-1].has_children = True
        self.stack.append(StackLevel(name, dict(attrs.items())))
        self.buffer = ""

    def characters(self, text):
        self.buffer += text

    def endElement(self, name):
        stack_id = list(i.name for i in self.stack)
        level = self.stack.pop()
        for s, out_name, f in f_map:
            if stack_match(s, stack_id):
                if out_name is None:
                    out_name = stack_id[-1]
                v = f(self.record_write, self.buffer, level.attrs, **level.data)
                if v is not None:
                    if isinstance(v, list):
                        if out_name in self.stack[-1].data:
                            self.stack[-1].data[out_name].extend(v)
                        else:
                            self.stack[-1].data[out_name] = v
                    elif isinstance(v, dict):
                        if out_name in self.stack[-1].data:
                            self.stack[-1].data[out_name] = dict(self.stack[-1].data, **v)
                        else:
                            self.stack[-1].data[out_name] = v
                    else:
                        self.stack[-1].data[out_name] = v
        # if not found:
        #     n = ",".join(stack_id)
        #     if n not in NOT_FOUND:
        #         NOT_FOUND.add(n)
        #         logging.warning("combiner for %s not found" % (",".join(stack_id)))
        self.buffer = ""


def parse_pubmed(handle, queue):
    handler = PubMedHandler(queue)
    parser = xml.sax.make_parser()
    parser.setContentHandler(handler)
    parser.parse(handle)


def name_clean(path):
    return os.path.basename(path).replace(".xml.gz", "")


def convert(path):
    out = []
    with gzip.open(path) as handle:
        parse_pubmed(handle, out)
    return out

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    parser = argparse.ArgumentParser()
    parser.add_argument("-l", action="store_true", default=False)
    parser.add_argument("-o", "--output", default="pubmed")
    parser.add_argument("--scan", default="source/pubmed/baseline")
    parser.add_argument("-N", default=1, type=int)

    parser.add_argument("files", nargs="*")
    args = parser.parse_args()

    if args.l:
        ftp = FTP('ftp.ncbi.nlm.nih.gov')
        ftp.login()
        for i in ftp.nlst("/pubmed/baseline/*.xml.gz"):
            name = os.path.basename(i)
            print(json.dumps({"url": "ftp://ftp.ncbi.nlm.nih.gov%s" % i, "name": name}))
        sys.exit(0)

    if len(args.files) == 0:
        args.files = glob(os.path.join(args.scan, "*.xml.gz"))

    with Pool(args.N) as pool:
        m = Manager()
        outputs = []
        for f in args.files:
            outputs.append( pool.apply_async(convert, (f,)) )

        emitter = JSONEmitter(directory=args.output, prefix="pubmed")
        for i in outputs:
            for o in i.get():
                p =  Publication(id=Publication.make_gid(o["url"]),
                                  url=o["url"], title=o["title"], abstract=o["abstract"],
                                  text="", date=o["date"], author=o["author"], citation=[],
                                  project_id=Project.make_gid("Reference"))
                emitter.emit_vertex(p)
        emitter.close()
