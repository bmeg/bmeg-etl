import csv
import sys
import os
from zipfile import ZipFile
import sqlite3
import json

# increase size of field
csv.field_size_limit(sys.maxsize)


def dict_factory(cursor, row):
    """sqlite utility"""
    d = {}
    for idx, col in enumerate(cursor.description):
        d[col[0]] = row[idx]
    return d


def load_table(conn, reader, table_name, call_back=None):
    """Loads sqlite table."""
    c = conn.cursor()
    # Create table
    c.execute(f'DROP TABLE IF EXISTS {table_name};')
    c.execute(f'CREATE TABLE {table_name} (id text, json text);')
    # print(f'CREATE TABLE {table_name} (id text, json text);')

    expected_count = 0
    for line in reader:
        id = list(line.values())[0]
        if call_back:
            line = call_back(line)
        serialized = json.dumps(line)
        # Insert a row of data
        c.execute(f"INSERT INTO {table_name}(id, json) VALUES (?, ?)", (id, serialized))
        expected_count += 1
    # Save (commit) the changes
    c.execute(f'CREATE INDEX {table_name}_idx on {table_name}(id);')
    conn.commit()
    c.execute(f'SELECT count(*) as "actual_count" from {table_name};')
    count = c.fetchone()
    count['expected_count'] = expected_count
    count['counts_ok'] = count['expected_count'] == count['actual_count']
    print(table_name, count)


def csv_field(line, fields):
    """Creates array from csv field, strips double quotes from fields in line."""
    for k in fields:
        line[k] = [a.replace('"', '').strip() for a in line[k].split(',') if a]
    return line


def csv_ids(line, fields):
    """Creates array from csv field, {name, id} where 'name (id)'"""
    for k in fields:
        # only one (or less ids in this string)
        if line[k].count('(PA') < 2:
            a = [line[k]]
        else:
            # start of string
            s = 0
            # fill this array with ids
            a = []
            # short hand
            d = line[k]
            while s < len(d):
                # find id
                e1 = d[s:].index('(PA')
                # find end of the id
                e2 = d[e1 + s:].index(')') + s
                e = e1 + e2 + 1
                # fill array with >name (PAXXXX)<
                a.append(d[s:e].replace('"', '').strip())
                # start at next string
                s = e + 1
        # fill this array with id objects
        _a = []
        for s in a:
            if s:
                _a.append({'name': before_parenthesis(s), 'id': between_parenthesis(s)})
        line[k] = _a
    return line


def csv_id(line, k):
    """Creates obj from csv field, {name, id} where 'name (id)'"""
    a = line[k].split('","')
    _a = []
    for s in a:
        if s:
            _a.append({'name': before_parenthesis(s), 'id': between_parenthesis(s)})
    line[k] = _a
    return line


def ontology_ids(line, fields):
    """Creates array from csv field, {name, id} where 'name (id)'"""
    for k in fields:
        a = line[k].split(',')
        _a = []
        for s in a:
            if s:
                _a.append({'name': between_parenthesis(s), 'id': before_parenthesis(s)})
        line[k] = _a
    return line


def load_all(conn, path):
    """Loads sqllite db."""
    with ZipFile(os.path.join(path, 'annotations.zip'), 'r') as zipObj:
        print('Extracting annotations')
        zipObj.extractall(path=path)
    with ZipFile(os.path.join(path, 'drugs.zip'), 'r') as zipObj:
        print('Extracting drugs')
        zipObj.extractall(path=path)
    with ZipFile(os.path.join(path, 'phenotypes.zip'), 'r') as zipObj:
        print('Extracting phenotypes')
        zipObj.extractall(path=path)
    with ZipFile(os.path.join(path, 'genes.zip'), 'r') as zipObj:
        print('Extracting genes')
        zipObj.extractall(path=path)
    with ZipFile(os.path.join(path, 'chemicals.zip'), 'r') as zipObj:
        print('Extracting chemicals')
        zipObj.extractall(path=path)

    #
    clinical_ann_metadata_keys = ['clinical_annotation_id', 'location', 'gene', 'level_of_evidence',
                                  'clinical_annotation_types', 'genotype_phenotype_ids', 'annotation_text',
                                  'variant_annotations_ids', 'variant_annotations', 'pubmed_ids', 'evidence_count',
                                  'related_chemicals', 'related_diseases', 'race', 'chromosome']
    fh = open(os.path.join(path, 'clinical_ann_metadata.tsv'), "r")
    tsv_in = csv.DictReader(fh, delimiter="\t", fieldnames=clinical_ann_metadata_keys, quoting=csv.QUOTE_NONE)
    next(tsv_in)  # skip 1st line

    #
    def clinical_ann_metadata_cb(line):  # callback before db write
        line = csv_field(line, ['genotype_phenotype_ids', 'variant_annotations_ids', 'pubmed_ids', 'clinical_annotation_types'])
        line = csv_id(line, 'gene')
        return csv_ids(line, ['related_diseases', 'related_chemicals'])
    load_table(conn, tsv_in, 'clinical_ann_metadata', clinical_ann_metadata_cb)
    fh.close()

    #
    clinical_ann_keys = ['genotype_phenotype_id', 'genotype', 'clinical_phenotype']
    fh = open(os.path.join(path, 'clinical_ann.tsv'), "r")
    tsv_in = csv.DictReader(fh, delimiter="\t", fieldnames=clinical_ann_keys, quoting=csv.QUOTE_NONE)
    next(tsv_in)  # skip 1st line
    load_table(conn, tsv_in, 'clinical_ann')
    fh.close()

    #
    var_drug_ann_keys = ['annotation_id', 'variant', 'gene', 'chemical', 'pubmed_id', 'phenotype_category', 'significance', 'notes', 'sentence', 'study_parameters', 'alleles', 'chromosome']
    fh = open(os.path.join(path, 'var_drug_ann.tsv'), "r")
    tsv_in = csv.DictReader(fh, delimiter="\t", fieldnames=var_drug_ann_keys, quoting=csv.QUOTE_NONE)
    next(tsv_in)  # skip 1st line

    def var_drug_ann_cb(line):  # callback before db write
        line = csv_field(line, ['pubmed_id', 'study_parameters'])
        return csv_ids(line, ['gene', 'chemical'])
    load_table(conn, tsv_in, 'var_drug_ann', var_drug_ann_cb)
    fh.close()

    #
    var_fa_ann_keys = ['annotation_id', 'variant', 'gene', 'chemical', 'pubmed_id', 'phenotype_category', 'significance', 'notes', 'sentence', 'study_parameters', 'alleles', 'chromosome']
    fh = open(os.path.join(path, 'var_fa_ann.tsv'), "r")
    tsv_in = csv.DictReader(fh, delimiter="\t", fieldnames=var_fa_ann_keys, quoting=csv.QUOTE_NONE)
    next(tsv_in)  # skip 1st line

    def var_fa_ann_cb(line):  # callback before db write
        line = csv_field(line, ['study_parameters'])
        return csv_ids(line, ['chemical', 'gene'])
    load_table(conn, tsv_in, 'var_fa_ann', var_fa_ann_cb)
    fh.close()

    #
    drugs_keys = ['pharm_gkb_accession_id', 'name', 'generic_names', 'trade_names', 'brand_mixtures', 'type',
                  'cross_references', 'smiles', 'in_chi', 'dosing_guideline', 'external_vocabulary',
                  'clinical_annotation_count', 'variant_annotation_count', 'pathway_count',
                  'vip_count', 'dosing_guideline_sources', 'top_clinical_annotation_level',
                  'top_fda_label_testing_level', 'top_any_drug_label_testing_level', 'label_has_dosing_info',
                  'has_rx_annotation', 'rx_norm_identifiers', 'atc_identifiers', 'pub_chem_compound_identifiers']
    fh = open(os.path.join(path, 'drugs.tsv'), "r")
    tsv_in = csv.DictReader(fh, delimiter="\t", fieldnames=drugs_keys, quoting=csv.QUOTE_NONE)
    next(tsv_in)  # skip 1st line

    def drugs_cb(line):  # callback before db write
        line = csv_field(line, ['rx_norm_identifiers', 'cross_references', 'trade_names', 'brand_mixtures'])
        return ontology_ids(line, ['external_vocabulary'])
    load_table(conn, tsv_in, 'drugs', drugs_cb)
    fh.close()

    #
    phenotypes_keys = ['pharm_gkb_accession_id', 'name', 'alternate_names', 'cross_references', 'external_vocabulary']
    fh = open(os.path.join(path, 'phenotypes.tsv'), "r")
    tsv_in = csv.DictReader(fh, delimiter="\t", fieldnames=phenotypes_keys, quoting=csv.QUOTE_NONE)
    next(tsv_in)  # skip 1st line

    def phenotypes_cb(line):  # callback before db write
        line = csv_field(line, ['alternate_names'])
        return ontology_ids(line, ['external_vocabulary', 'cross_references'])
    load_table(conn, tsv_in, 'phenotypes', phenotypes_cb)
    fh.close()

    #
    genes_keys = ['pharm_gkb_accession_id', 'ncbi_gene_id', 'hgnc_id', 'ensembl_id', 'name', 'symbol',
                  'alternate_names', 'alternate_symbols', 'is_vip', 'has_variant_annotation', 'cross_references',
                  'has_cpic_dosing_guideline', 'chromosome',
                  'chromosomal_start_grch37_p13', 'chromosomal_stop_grch37_p13', 'chromosomal_start_grch38_p7', 'chromosomal_stop_grch38_p7']
    fh = open(os.path.join(path, 'genes.tsv'), "r")
    tsv_in = csv.DictReader(fh, delimiter="\t", fieldnames=genes_keys, quoting=csv.QUOTE_NONE)
    next(tsv_in)  # skip 1st line

    def genes_cb(line):  # callback before db write
        line = csv_field(line, ['alternate_names', 'cross_references', 'alternate_symbols'])
        return line
    load_table(conn, tsv_in, 'genes', genes_cb)
    fh.close()

    #
    var_pheno_ann_keys = ['annotation_id', 'variant', 'gene', 'chemical', 'pubmed_id', 'phenotype_category', 'significance', 'notes', 'sentence', 'study_parameters', 'alleles', 'chromosome']
    fh = open(os.path.join(path, 'var_pheno_ann.tsv'), "r")
    tsv_in = csv.DictReader(fh, delimiter="\t", fieldnames=var_pheno_ann_keys, quoting=csv.QUOTE_NONE)
    next(tsv_in)  # skip 1st line

    def var_pheno_ann_cb(line):  # callback before db write
        line = csv_ids(line, ['gene', 'chemical'])
        return line
    load_table(conn, tsv_in, 'var_pheno_ann', var_pheno_ann_cb)
    fh.close()

    #
    chemicals_keys = ['pharm_gkb_accession_id', 'name', 'generic_names', 'trade_names', 'brand_mixtures', 'type', 'cross_references',
                      'smiles', 'in_chi', 'dosing_guideline', 'external_vocabulary', 'clinical_annotation_count',
                      'variant_annotation_count', 'pathway_count', 'vip_count', 'dosing_guideline_sources', 'top_clinical_annotation_level',
                      'top_fda_label_testing_level', 'top_any_drug_label_testing_level', 'label_has_dosing_info',
                      'has_rx_annotation', 'rx_norm_identifiers', 'atc_identifiers', 'pub_chem_compound_identifiers']
    fh = open(os.path.join(path, 'chemicals.tsv'), "r")
    tsv_in = csv.DictReader(fh, delimiter="\t", fieldnames=chemicals_keys, quoting=csv.QUOTE_NONE)
    next(tsv_in)  # skip 1st line

    def chemicals_cb(line):  # callback before db write
        line = csv_field(line, ['cross_references', 'pub_chem_compound_identifiers', 'rx_norm_identifiers', 'generic_names', 'trade_names', 'brand_mixtures'])
        return ontology_ids(line, ['external_vocabulary'])
    load_table(conn, tsv_in, 'chemicals', chemicals_cb)
    fh.close()


def between_parenthesis(s):
    """Parses string, returns XXX in '.*(XXX).*', Strips blanks. )"""
    if '(' not in s:
        return None
    return s[s.rfind("(") + 1: s.rfind(")")].strip()


def before_parenthesis(s):
    """Parses string, returns XXX in 'XXX (.*', Strips blanks. )"""
    if '(' not in s:
        return s
    s = s[0:s.rfind("(")].strip()
    if s[0] == ',':
        s = s[1:]
    s = s.replace('"', '')
    return s


def transform(source_path="source/pharmgkb",
              emitter_prefix=None,
              emitter_directory='pharmgkb'):
    """"""
    path = os.path.abspath(source_path)
    db_file = os.path.join(path, 'work.sqlite')
    if os.path.exists(db_file):
        print(f'{db_file} exists, rebuilding')
        os.remove(db_file)
    conn = sqlite3.connect(db_file)
    conn.row_factory = dict_factory
    load_all(conn, path)

    conn.close()


if __name__ == '__main__':  # pragma: no cover
    transform()
