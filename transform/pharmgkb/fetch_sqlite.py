import os
import sqlite3
import json


def dict_factory(cursor, row):
    """sqlite utility"""
    d = {}
    for idx, col in enumerate(cursor.description):
        d[col[0]] = row[idx]
    return d


def load_ids(conn, line, id_field, table_name, alternate_table_names=None):
    """Hydrates id_field csv from table_name."""
    _a = []
    c = conn.cursor()
    for id in line[id_field]:
        v = c.execute(f'select * from {table_name} where id = ? ;', (id['id'],)).fetchone()
        if not v and alternate_table_names:
            for alternate_table_name in alternate_table_names:
                v = c.execute(f'select * from {alternate_table_name} where id = ? ;', (id['id'],)).fetchone()
                if v:
                    break
        if v:
            _a.append(json.loads(v['json']))
        else:
            print(f'WARNING >{id}< not found in [{table_name}, {alternate_table_names}] clinical_ann_metadata.id {list(line.values())[0]}')
    line[id_field] = _a
    return line


def load_ids_scalar(conn, line, id_field, table_name, alternate_table_names=None):
    """Hydrates id csv from table_name."""
    _a = []
    c = conn.cursor()
    for id in line[id_field]:
        v = c.execute(f'select * from {table_name} where id = ? ;', (id,)).fetchone()
        if not v and alternate_table_names:
            for alternate_table_name in alternate_table_names:
                v = c.execute(f'select * from {alternate_table_name} where id = ? ;', (id,)).fetchone()
                if v:
                    break
        if v:
            _a.append(json.loads(v['json']))
        else:
            print(f'WARNING >{id}< not found in [{table_name}, {alternate_table_names}] clinical_ann_metadata.id {list(line.values())[0]}')
    line[id_field] = _a
    return line


def load_id(conn, id, table_name, alternate_table_names=None):
    """Load a single record from table_name."""
    c = conn.cursor()
    v = c.execute(f'select * from {table_name} where id = ? ;', (id,)).fetchone()
    if not v and alternate_table_names:
        for alternate_table_name in alternate_table_names:
            v = c.execute(f'select * from {alternate_table_name} where id = ? ;', (id,)).fetchone()
            if v:
                break
    # if not v:
    #     print(f'WARNING >{id}< not found in [{table_name}, {alternate_table_names}]')
    return v


def sqlite_connection(source_path):
    """Create connection to sqlite."""
    path = os.path.abspath(source_path)
    db_file = os.path.join(path, 'work.sqlite')
    assert os.path.exists(db_file), f'{db_file} missing'
    conn = sqlite3.connect(db_file)
    conn.row_factory = dict_factory
    return conn


def features(annotation, conn, include_haplotypes=False):
    """Returns all snp ids in annotation, drills down through haplotypes."""
    for variant_annotation in annotation['variant_annotations_ids']:
        variant = variant_annotation['variant']
        if variant.startswith('rs'):
            yield variant
        else:
            # haplotypes have messy representation CSV with embedded commas :-(
            # use havested table to verify
            miss = True
            for l in [l.strip() for l in variant.split(',')]:
                obj = load_id(conn, l, 'haplotypes', ['variants'])
                if obj:
                    miss = False
                    if include_haplotypes:
                        yield l
                    d = json.loads(obj['json'])
                    if 'variants' in d:
                        # is haplotype
                        for v in [v for v in d['variants'] if v.startswith('rs')]:
                            yield v
            if miss and ',' in variant:
                obj = load_id(conn, variant.strip(), 'haplotypes', ['variants'])
                if obj:
                    miss = False
                    d = json.loads(obj['json'])
                    if include_haplotypes:
                        yield variant.strip()
                    if 'variants' in d:
                        # is haplotype
                        for v in [v for v in d['variants'] if v.startswith('rs')]:
                            yield v
            if miss:
                print(f"WARNING >{variant}< not found")


def fetch_feature(conn, id):
    """Returns genomic coordinates."""
    o = load_id(conn, id, 'refsnp')
    if not o:
        print(f"WARNING {id} not found?")
        return None
    o = json.loads(o['json'])
    if 'mappings' not in o:
        print(f"WARNING {id} has no mappings? {o}")
        return None
    hg37 = [g for g in o.get('mappings', []) if g['assembly_name'] == 'GRCh37']
    if len(hg37) == 0:
        print(f"WARNING {id} has not hg37? {o}")
        return None
    hg37[0]['rsid'] = id
    gene_symbols = []
    if 'phenotypes' in o:
        for p in o['phenotypes']:
            if p['genes']:
                gene_symbols.extend(p['genes'].split(','))
    hg37[0]['gene_symbols'] = list(set(gene_symbols))
    return hg37[0]


def fetch_all(source_path=None, conn=None, include_features=True):
    """Fetches all clinical annotations and child nodes from db, yields"""
    close_after = True
    if conn:
        c = conn.cursor()
        close_after = False
    else:
        conn = sqlite_connection(source_path)
        c = conn.cursor()
    for clinical_ann_metadata in c.execute('select * from clinical_ann_metadata').fetchall():
        clinical_ann_metadata = json.loads(clinical_ann_metadata['json'])
        clinical_ann_metadata = load_ids(conn, clinical_ann_metadata, 'related_chemicals', 'drugs', ['chemicals'])
        clinical_ann_metadata = load_ids(conn, clinical_ann_metadata, 'related_diseases', 'phenotypes')
        clinical_ann_metadata = load_ids(conn, clinical_ann_metadata, 'gene', 'genes')
        clinical_ann_metadata = load_ids_scalar(conn, clinical_ann_metadata, 'variant_annotations_ids', 'var_drug_ann', ['var_pheno_ann', 'var_fa_ann'])
        if include_features:
            clinical_ann_metadata['feature_names'] = [f for f in set(features(clinical_ann_metadata, conn, include_haplotypes=True))]
            clinical_ann_metadata['features'] = [fetch_feature(conn, f) for f in set(features(clinical_ann_metadata, conn))]
            clinical_ann_metadata['features'] = [f for f in clinical_ann_metadata['features'] if f]  # clean up no finds
        yield(clinical_ann_metadata)
    if close_after:
        conn.close()
