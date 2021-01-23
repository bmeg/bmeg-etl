#!/usr/bin/env python

import os
import sys
import sqlite3
import argparse


def chembl_synonym_search(conns, name):
    c = conns.chembl.cursor()
    c.execute("select * from MOLECULE_SYNONYMS where UPPER(synonyms) == ?", (name.upper(),))
    names = list(i[0] for i in c.description)
    ids = set()
    for row in c:
        data = dict(zip(names,row))
        i = data['molregno']
        ids.add(i)

    out = []
    for i in ids:
        c.execute("select chembl_id from MOLECULE_DICTIONARY where MOLREGNO=?",(i,))
        names = list(i[0] for i in c.description)
        for row in c:
            data = dict(zip(names,row))
            out.append(data["chembl_id"])
    return out


def chembl_search_name(conns, name):
    c = conns.chembl.cursor()
    out = []
    c.execute("select chembl_id from MOLECULE_DICTIONARY where UPPER(pref_name)=?",[name.upper()])
    for row in c:
        out.append(row[0])
    return out


def get_parent(conns, id):
    c = conns.chembl.cursor()
    molregno = ""
    for row in c.execute("select molregno from MOLECULE_DICTIONARY where chembl_id=?", [id]):
        molregno = row[0]

    parent_molregno = ""
    for row in c.execute("select parent_molregno from molecule_hierarchy where molregno=?", [molregno]):
        parent_molregno = row[0]

    parent_id = None
    for row in c.execute("select chembl_id from molecule_dictionary where molregno=?", [parent_molregno]):
        parent_id = row[0]

    return parent_id

def search_pubchem2chembl(conns, id):
    c = conns.pubchem.cursor()
    out = []
    for row in c.execute('SELECT name FROM synonym where id = ? and name like "CHEMBL%"', [id.replace("CID", "")]):
        out.append(row[0])
    return out

def search_pubchem_alias(conns, term):
    c = conns.pubchem.cursor()
    out = []
    for row in c.execute('SELECT id FROM synonym where name = ?', [term.upper()]):
        out.append("CID%s" % (row[0]))
    return out

def pubchem_synonym_lookup(conns, args):
    c = conns.pubchem.cursor()
    with open(args.input) as handle:
        for line in handle:
            n=line.rstrip()
            if not n.startswith("CHEMBL"):
                sn = n.upper()
                if args.rd:
                    sn = sn.replace("-", " ")
                for row in c.execute('SELECT name,id FROM synonym where name = ?', [sn]):
                    print("%s\tCID%s" % (n,row[1]))

def pubchem2chembl_lookup(conns, args):
    c = conns.pubchem.cursor()
    with open(args.file) as handle:
        if args.mapfile:
            for s in handle:
                term, name = s.rstrip().split("\t")
                if name.startswith("CID"):
                    for row in c.execute('SELECT id,name FROM synonym where id = ? and name like "CHEMBL%"', [name.replace("CID", "")]):
                        print("%s\t%s" % (term,row[1]))
        else:
            for s in handle:
                name = s.rstrip()
                if name.startswith("CID"):
                    for row in c.execute('SELECT id,name FROM synonym where id = ? and name like "CHEMBL%"', [name.replace("CID", "")]):
                        print("CID%s\t%s" % (row[0],row[1]))


def chembl_name_lookup(conns, args):
    names = set()
    with open(args.input) as handle:
        for line in handle:
            n = line.rstrip()
            names.add(n)

    for name in names:
        if args.rd:
            o = search_name(conns, name.replace("-", " "))
        else:
            o = search_name(conns, name)
        for i in o:
            print("%s\t%s" % (name, i))



def chembl_synonym_lookup(conns, args):
    names = set()
    with open(args.input) as handle:
        for line in handle:
            n = line.rstrip()
            names.add(n)

    for name in names:
        if args.rd:
            o = chembl_synonym_search(conns, name.replace("-", " "))
        else:
            o = chembl_synonym_search(conns, name)
        for i in o:
            print("%s\t%s" % (name, i))


def dedup(conns, args):
    termMap = {}
    with open(args.input) as handle:
        for line in handle:
            term, id = line.rstrip().split("\t")
            termMap[term] = termMap.get(term, set()).union( set([id]) )

    for term, idSet in termMap.items():
        if len(idSet) > 1:
            parentSet = set()
            for i in idSet:
                p = get_parent(conns, i)
                if p:
                    parentSet.add(p)
            for i in parentSet:
                print("%s\t%s" % (term, i))
        else:
            print("%s\t%s" % (term, list(idSet)[0]))


def search(conns, args):
    nameMap = conns.get_name_table()
    with open(args.input) as handle:
        for line in handle:
            term = line.rstrip()
            found = False

            if term.startswith("CHEMBL"):
                print("%s\t%s" % (term,term) )
                found = True
            if not found:
                if term.startswith("CID"):
                    for i in search_pubchem2chembl(conns, term):
                        print("%s\t%s" % (term,i) )
                        found = True
            if not found:
                if term.upper() in nameMap:
                    for i in nameMap[term.upper()]:
                        print("%s\t%s" % (term,i) )
                        found = True
            if not found:
                t = term.upper().replace("-", " ")
                if t in nameMap:
                    for i in nameMap[t]:
                        print("%s\t%s" % (term,i) )
                        found = True
            if not found:
                for i in chembl_synonym_search(conns, term):
                    print("%s\t%s" % (term,i) )
                    found = True
            if not found:
                for i in chembl_synonym_search(conns, term.replace("-", " ")):
                    print("%s\t%s" % (term,i) )
                    found = True

            if not found:
                for cid in search_pubchem_alias(conns, term):
                    for i in search_pubchem2chembl(conns, cid):
                        print("%s\t%s" % (term,i) )
                        found = True

            if not found:
                for cid in search_pubchem_alias(conns, term.replace("-", " ")):
                    for i in search_pubchem2chembl(conns, cid):
                        print("%s\t%s" % (term,i) )
                        found = True


class Connections:
    def __init__(self, chembl_conn, pubchem_conn):
        self.chembl = chembl_conn
        self.pubchem = pubchem_conn

    def get_name_table(self):
        c = self.chembl.cursor()
        out = {}
        for row in c.execute("SELECT pref_name,chembl_id from MOLECULE_DICTIONARY"):
            if row[0] is not None:
                name = row[0].upper()
                out[name] = out.get(name, set()).union([row[1]])
        return out

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--database", help="Database File")
    parser.add_argument("-s", "--synonym", help="Synonym File")

    subparsers = parser.add_subparsers(help='sub-command help')

    chembl_synonym_parser = subparsers.add_parser('chembl-synonym')
    chembl_synonym_parser.add_argument('input')
    chembl_synonym_parser.add_argument("--rd", action="store_true", default=False, help="Remove Dashes")
    chembl_synonym_parser.set_defaults(func=chembl_synonym_lookup)

    chembl_name_parser = subparsers.add_parser('chembl-name')
    chembl_name_parser.add_argument('input')
    chembl_name_parser.add_argument("--rd", action="store_true", default=False, help="Remove Dashes")
    chembl_name_parser.set_defaults(func=chembl_name_lookup)

    dedup_parser = subparsers.add_parser('dedup')
    dedup_parser.add_argument('input')
    dedup_parser.set_defaults(func=dedup)

    search_parser = subparsers.add_parser('search')
    search_parser.add_argument('input')
    search_parser.set_defaults(func=search)

    pubchem_synonym_parser = subparsers.add_parser('pubchem-synonym')
    pubchem_synonym_parser.add_argument("--rd", action="store_true", default=False, help="Remove Dashes")
    pubchem_synonym_parser.add_argument('input')
    pubchem_synonym_parser.set_defaults(func=pubchem_synonym_lookup)

    pubchem2chembl_parser = subparsers.add_parser('pubchem2chembl')
    pubchem2chembl_parser.add_argument('input')
    pubchem2chembl_parser.add_argument('-m', '--mapfile', action="store_true", default=False)
    pubchem2chembl_parser.set_defaults(func=pubchem2chembl_lookup)

    args = parser.parse_args()

    chembl_conn = sqlite3.connect(args.database)

    if not os.path.exists(args.synonym + ".sqlite"):
        pubchem_conn = sqlite3.connect(args.synonym + ".sqlite")
        c = pubchem_conn.cursor()
        c.execute("CREATE TABLE synonym(id int, name text)")
        c.execute("CREATE INDEX name_idx on synonym(name)")
        c.execute("CREATE INDEX id_idx on synonym(id)")
        with gzip.open(args.synonym, "rt") as handle:
            for line in handle:
                id, term = line.rstrip().split("\t")
                if term.endswith(";"):
                    term = term[:-1]
                c.execute("INSERT INTO synonym VALUES(?,?)", [id,term.upper()] )
        pubchem_conn.commit()
    else:
        pubchem_conn = sqlite3.connect(args.synonym + ".sqlite")

    conns = Connections(chembl_conn, pubchem_conn)

    args.func(conns, args)
