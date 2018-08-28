#!/usr/bin/env python3

import argparse
import os
import tarfile
import requests

from xml.dom.minidom import parseString

from bmeg.vertex import PFAMFamily, PFAMClan, GeneOntologyTerm
from bmeg.edge import GeneOntologyAnnotation, PFAMClanMember
from bmeg.emitter import JSONEmitter


def getText(nodelist):
    rc = []
    for node in nodelist:
        if node.nodeType == node.TEXT_NODE:
            rc.append(node.data)
        elif node.nodeType == node.CDATA_SECTION_NODE:
            rc.append(node.data)
    return ''.join(rc)


def dom_scan(node, query):
    stack = query.split("/")
    if node.localName == stack[0] or stack[0] == "*":
        return dom_scan_iter(node, stack[1:], [stack[0]])


def dom_scan_iter(node, stack, prefix):
    if len(stack):
        for child in node.childNodes:
            if child.nodeType == child.ELEMENT_NODE:
                if child.localName == stack[0]:
                    for out in dom_scan_iter(child, stack[1:],
                                             prefix + [stack[0]]):
                        yield out
                elif '*' == stack[0]:
                    for out in dom_scan_iter(child, stack[1:],
                                             prefix + [child.localName]):
                        yield out
    else:
        if node.nodeType == node.ELEMENT_NODE:
            yield node, prefix, dict(node.attributes.items()), getText(
                node.childNodes)
        elif node.nodeType == node.TEXT_NODE:
            yield node, prefix, None, getText(node.childNodes)


def xml_transform(dom, emit):
    for elem, path, attr, _ in dom_scan(dom, "*/pfam/entry"):
        pfam_id = attr["id"]
        pfam_acc = attr["accession"]
        pfam_type = attr["entry_type"]
        comments = ""
        for _, _, _, txt in dom_scan(elem, "*/comment"):
            comments = txt
        description = ""
        for _, _, _, txt in dom_scan(elem, "*/description"):
            description = txt
        go_terms = []
        for _, _, attr, _ in dom_scan(elem, "*/go_terms/category/term"):
            go_terms.append(attr["go_id"])

        clans = []
        for _, _, attr, _ in dom_scan(elem, "*/clan_membership"):
            clans.append(attr['clan_acc'])
        """
        out.id = pfam_id
        out.accession = pfam_acc
        out.type = pfam_type
        out.go_terms.extend(go_terms)
        out.clans.extend(clans)
        out.description = description.strip()
        out.comments = comments.strip()
        """

        out = PFAMFamily(
            pfam_id=pfam_id, accession=pfam_acc, type=pfam_type,
            description=description.strip(), comments=comments.strip()
        )
        emit.emit_vertex(out)
        for g in go_terms:
            emit.emit_edge(
                GeneOntologyAnnotation(evidence="NA", title="", references=[]),
                from_gid=GeneOntologyTerm.make_gid(g),
                to_gid=out.gid()
            )
        for c in clans:
            emit.emit_edge(
                PFAMClanMember(),
                from_gid=PFAMClan.make_gid(c),
                to_gid=out.gid()
            )


def run_file(args):
    emitter = JSONEmitter(args.output_prefix)

    if args.in_archive is not None:
        tar = tarfile.open(args.in_archive, "r:gz")
        for member in tar.getmembers():
            handle = tar.extractfile(member)
            if handle is not None:
                dom = parseString(handle.read())
                xml_transform(dom, emitter)

    for f in args.files:
        with open(f) as handle:
            dom = parseString(handle.read())
            xml_transform(dom, emitter)


def get_acc_list():
    out = []
    handle = requests.get("http://pfam.xfam.org/families?output=text")
    for line in handle.iter_lines():
        line = line.decode()
        if not line.startswith("#"):
            row = line.split("\t")
            if len(row[0]) > 1:
                out.append(row[0])
    return out


def run_list(args):
    for i in get_acc_list():
        print(i)


def run_download(args):
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    ids = args.ids
    if args.all:
        ids = get_acc_list()

    if args.archive:
        tar = tarfile.open(
            os.path.join(args.output_dir, "pfam.tar.gz"), "w:gz")

    for i in ids:
        print(i)
        handle = requests.get(
            "http://pfam.xfam.org/family?output=xml&acc=%s" % i)
        txt = handle.text
        f = os.path.join(args.output_dir, "%s.xml" % (i))
        with open(f, "w") as handle:
            handle.write(txt)

        if args.archive:
            tar.add(f, arcname=os.path.basename(f))

    if args.archive:
        tar.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    subparser = parser.add_subparsers()

    parser_list = subparser.add_parser("list")
    parser_list.set_defaults(func=run_list)

    parser_download = subparser.add_parser("download")
    parser_download.set_defaults(func=run_download)
    parser_download.add_argument(
        "--all", "-a", action="store_true", default=False)
    parser_download.add_argument(
        "--output-dir",
        "-o",
        default="./",
        help="diretory in which  to create output files")
    parser_download.add_argument(
        "--archive",
        "-A",
        action="store_true",
        default=False,
        help="create a gzipped TAR archive of all the output files")
    parser_download.add_argument("ids", nargs="*")

    parser_file = subparser.add_parser("transform")
    parser_file.set_defaults(func=run_file)
    parser_file.add_argument(
        "--in-archive", help="gzipped tar archive containing PFAM XML files")
    parser_file.add_argument(
        "--output-prefix", "-o", default="pfam_output", help="output file")
    parser_file.add_argument("files", nargs="*", help="PFAM XML file(s)")

    args = parser.parse_args()
    args.func(args)
