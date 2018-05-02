#!/usr/bin/env python

# parses file found at ftp://ftp.ebi.ac.uk/pub/databases/genenames/hgnc_complete_set.txt.gz

# example usage -------------------

# python -m convert.hugo.convert-hugo --hugo ~/Data/hugo/hgnc_complete_set.txt --out ~/Data/hugo/hugo.json

# example hugo entry --------------------

# {'Accession Numbers': 'AF271790',
#  'Approved Name': 'APOBEC1 complementation factor',
#  'Approved Symbol': 'A1CF',
#  'CCDS IDs': 'CCDS7241, CCDS7242, CCDS7243, CCDS73133',
#  'Chromosome': '10q21.1',
#  'Date Approved': '2007-11-23',
#  'Date Modified': '2014-11-19',
#  'Date Name Changed': '',
#  'Date Symbol Changed': '',
#  'Ensembl Gene ID': 'ENSG00000148584',
#  'Ensembl ID (supplied by Ensembl)': 'ENSG00000148584',
#  'Entrez Gene ID': '29974',
#  'Entrez Gene ID (supplied by NCBI)': '29974',
#  'Enzyme IDs': '',
#  'Gene Family Tag': 'RBM',
#  'Gene family description': 'RNA binding motif (RRM) containing',
#  'HGNC ID': 'HGNC:24086',
#  'Locus Group': 'protein-coding gene',
#  'Locus Specific Databases': '',
#  'Locus Type': 'gene with protein product',
#  'Mouse Genome Database ID': 'MGI:1917115',
#  'Mouse Genome Database ID (supplied by MGI)': 'MGI:1917115',
#  'Name Synonyms': '',
#  'OMIM ID (supplied by NCBI)': '',
#  'Previous Names': '',
#  'Previous Symbols': '',
#  'Primary IDs': '',
#  'Pubmed IDs': '11815617, 11072063',
#  'Rat Genome Database ID (supplied by RGD)': 'RGD:619834',
#  'Record Type': 'Standard',
#  'RefSeq (supplied by NCBI)': 'XM_011539729',
#  'RefSeq IDs': 'NM_014576',
#  'Secondary IDs': '',
#  'Specialist Database IDs': ', , , , , , , , , , A1CF, , , , , , , ',
#  'Specialist Database Links': '<!--,--> <!--,--> <!--,--> <!--,--> <!--,--> <!--,--> <!--,--> <!--,--> <!--,--> <!--,--> <a href="http://cancer.sanger.ac.uk/cosmic/gene/overview?ln=A1CF">COSMIC</a><!--,--> <!--,--> <!--,--> <!--,--> <!--,--> <!--,--> <!--,--> ',
#  'Status': 'Approved',
#  'Synonyms': 'ACF, ASP, ACF64, ACF65, APOBEC1CF',
#  'UCSC ID (supplied by UCSC)': 'uc001jjj.3',
#  'UniProt ID (supplied by UniProt)': 'Q9NQ94',
#  'VEGA IDs': 'OTTHUMG00000018240',
#  'Vega ID (supplied by Vega)': 'OTTHUMG00000018240'}


import sys
import csv
import re
from pprint import pprint
import argparse
import json
from google.protobuf import json_format

from bmeg import genome_pb2, nlp_pb2, record

comma_match = r'\s*,\s*'

def split_zip(r, line, keys):
    if len(line[keys[0]]):
        return zip(*[re.split(r, line[key]) for key in keys])
    else:
        return []
    
class GeneGenerator(record.RecordGenerator):
    def __init__(self, family_gid, pubmed_gid):
        super(GeneGenerator, self).__init__('Gene')
        self.family_gid = family_gid
        self.pubmed_gid = pubmed_gid

    def schema(self):
        return genome_pb2.Gene()

    def gid(self, data):
        return 'gene:' + data['Approved Symbol']

    def update(self, gene, data):
        gene.symbol = data['Approved Symbol']
        gene.description = data['Approved Name']
        gene.chromosome = data['Chromosome']
        gene.accession = data['Accession Numbers']
        gene.refseq = data['RefSeq IDs']

        for family in data['families']:
            record.append_unique(gene.inFamily, self.family_gid(family))

        for citation in data['citations']:
            record.append_unique(gene.citedFrom, self.pubmed_gid(citation))

        return gene

class GeneSynonymGenerator(record.RecordGenerator):
    def __init__(self, key, gene_gid, database_gid):
        super(GeneSynonymGenerator, self).__init__('GeneSynonym')
        self.key = key
        self.gene_gid = gene_gid
        self.database_gid = database_gid

    def schema(self):
        return genome_pb2.GeneSynonym()

    def gid(self, data):
        return 'geneSynonym:' + data[self.key]

    def update(self, synonym, data):
        synonym.symbol = data[self.key]
        record.append_unique(synonym.synonymFor, self.gene_gid(data))
        record.append_unique(synonym.inDatabase, self.database_gid(data.get('database') or 'hugo'))
        return synonym

class GeneDatabaseGenerator(record.RecordGenerator):
    def __init__(self):
        super(GeneDatabaseGenerator, self).__init__('GeneDatabase')

    def schema(self):
        return genome_pb2.GeneDatabase()

    def gid(self, data):
        return 'geneDatabase:' + data

    def update(self, database, data):
        database.name = data
        return database

class GeneFamilyGenerator(record.RecordGenerator):
    def __init__(self):
        super(GeneFamilyGenerator, self).__init__('GeneFamily')

    def schema(self):
        return genome_pb2.GeneFamily()

    def gid(self, data):
        return 'geneFamily:' + data[0]

    def update(self, family, data):
        family.tag = data[0]
        family.description = data[1]
        return family

class PubmedGenerator(record.RecordGenerator):
    def __init__(self):
        super(PubmedGenerator, self).__init__('Pubmed')
        
    def schema(self):
        return nlp_pb2.Pubmed()

    def gid(self, data):
        return 'pubmed:' + data

    def update(self, pubmed, data):
        pubmed.pmid = data
        return pubmed

def initial_state():
    pubmed_generator = PubmedGenerator()
    database_generator = GeneDatabaseGenerator()
    family_generator = GeneFamilyGenerator()
    gene_generator = GeneGenerator(family_generator.gid, pubmed_generator.gid)

    types = ['Gene', 'GeneSynonym', 'GeneDatabase', 'GeneFamily', 'Pubmed']

    synonyms = [
        ('hugo', 'Synonyms'),
        ('ensembl', 'Ensembl Gene ID'),
        ('entrez', 'Entrez Gene ID'),
        ('hgnc', 'HGNC ID'),
        ('mgd', 'Mouse Genome Database ID'),
        ('rgd', 'Rat Genome Database ID (supplied by RGD)'),
        ('refseq', 'RefSeq (supplied by NCBI)'),
        ('ucsc', 'UCSC ID (supplied by UCSC)'),
        ('uniprot', 'UniProt ID (supplied by UniProt)'),
        ('vega', 'VEGA IDs')
    ]

    state = {
        'types': types,
        'synonyms': synonyms,
        'generators': {
            'gene': gene_generator,
            'family': family_generator,
            'database': database_generator,
            'pubmed': pubmed_generator
        }
    }

    for type in types:
        state[type] = {}

    state['generators']['identity'] = GeneSynonymGenerator('Approved Symbol', gene_generator.gid, database_generator.gid)
    for synonym in synonyms:
        state['generators'][synonym[0]] = GeneSynonymGenerator(synonym[1], gene_generator.gid, database_generator.gid)
        state['generators']['database'].find(state, synonym[0])

    return state

def split_append(state, line, key, transform=lambda l, p: p):
    records = []
    if len(line[key]):
        for part in re.split(comma_match, line[key]):
            record = transform(line, part)
            records.append(record)

    return records

def synonym_append(state, line, key, generator):
    def trigger_generator(line, part):
        return state['generators'][generator].find(state, {
            key: part,
            'database': generator,
            'Approved Symbol': line['Approved Symbol']
        })

    return split_append(state, line, key, trigger_generator)

def convert_line(state, line):
    if line['Status'] == 'Approved':
        symbol = line['Approved Symbol']
        print(symbol)

        line['families'] = split_zip(comma_match, line, ['Gene Family Tag', 'Gene family description'])
        line['citations'] = split_append(state, line, 'Pubmed IDs')
        line['synonyms'] = split_append(state, line, 'Synonyms')

        gene = state['generators']['gene'].find(state, line)

        state['generators']['identity'].find(state, {'database': 'hugo', 'Approved Symbol': symbol})
        for synonym in state['synonyms']:
            synonym_append(state, line, synonym[1], synonym[0])

        for family in line['families']:
            state['generators']['family'].find(state, family)

        for citation in line['citations']:
            state['generators']['pubmed'].find(state, citation)

    return state

def message_to_json(message):
    json = json_format.MessageToJson(message)
    return re.sub(r' +', ' ', json.replace('\n', ''))

def output_state(state, prefix):
    for type in state['types']:
        json = []
        for key in state[type]:
            message = message_to_json(state[type][key])
            json.append(message)

        out = '\n'.join(json)
        outhandle = open(prefix + '.' + type + '.json', 'w')
        outhandle.write(out)
        outhandle.close()

def convert_hugo(file):
    state = initial_state()
    reader = csv.DictReader(file, delimiter='\t')
    for line in reader:
        state = convert_line(state, line)

    return state

def parse_args(args):
    args = args[1:]
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('--hugo', type=str, help='path to the hugo source file')
    parser.add_argument('--out', type=str, default='hugo', help='Path to output files (.json or .pbf_ext)')
    parser.add_argument('--format', type=str, default='json', help='Format of output: json or pbf (binary)')

    return parser.parse_args(args)

if __name__ == '__main__':
    options = parse_args(sys.argv)
    state = convert_hugo(open(options.hugo))
    output_state(state, options.out)