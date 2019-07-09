#!/usr/bin/env python3

from glob import glob

from xml.dom.minidom import parseString

from bmeg import (PfamFamily, PfamClan, GeneOntologyTerm, Project,
                  GeneOntologyTerm_PfamFamilies_PfamFamily,
                  PfamClan_PfamFamilies_PfamFamily)

from bmeg.emitter import JSONEmitter
from bmeg.ioutils import read_tsv


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

        out = PfamFamily(
            id=PfamFamily.make_gid(pfam_acc),
            pfam_id=pfam_id,
            accession=pfam_acc,
            type=pfam_type,
            description=description.strip(),
            comments=comments.strip(),
            project_id=Project.make_gid("Reference")
        )
        emit.emit_vertex(out)
        for g in go_terms:
            emit.emit_edge(
                GeneOntologyTerm_PfamFamilies_PfamFamily(
                    from_gid=GeneOntologyTerm.make_gid(g),
                    to_gid=out.gid(),
                    data={'evidence': None, 'title': None, 'references': None}
                ),
                emit_backref=True
            )

        for c in clans:
            emit.emit_edge(
                PfamClan_PfamFamilies_PfamFamily(
                    from_gid=PfamClan.make_gid(c),
                    to_gid=out.gid()
                ),
                emit_backref=True
            )


def transform(pfam_xmls="source/pfam/xmls/*.xml",
              clans_file="source/pfam/clans.tsv",
              emitter_prefix=None,
              emitter_directory='pfam'):

    emitter = JSONEmitter(directory=emitter_directory, prefix=emitter_prefix)

    for f in glob(pfam_xmls):
        with open(f) as handle:
            dom = parseString(handle.read())
            xml_transform(dom, emitter)

    tsv_in = read_tsv(clans_file, fieldnames=["accession", "id", "description"])
    for line in tsv_in:
        # accession	id	description
        c = PfamClan(
            id=PfamClan.make_gid(line["accession"]),
            accession=line["accession"],
            clan_id=line["id"],
            description=line["description"],
            project_id=Project.make_gid("Reference")
        )
        emitter.emit_vertex(c)

    emitter.close()


if __name__ == '__main__':  # pragma: no cover
    transform()
