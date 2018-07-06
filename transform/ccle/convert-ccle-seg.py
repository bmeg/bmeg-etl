#!/usr/bin/env python

# curl -o CCLE_copynumber_2013-12-03.seg.txt "https://portals.broadinstitute.org/ccle/downloadFile/DefaultSystemRoot/exp_10/ds_20/CCLE_copynumber_2013-12-03.seg.txt?downloadff=true&fileId=17597"

import re
import sys
import gzip
import json
import bmeg.cna_pb2 as cna
from google.protobuf import json_format

from bx.intervals.intersection import Intersecter, Interval

reAttr = re.compile(r'([^ ]+) \"(.*)\"')

class GTFLine:
    def __init__(self, seqname, source, feature, start, end, score, strand, frame, attr):
        self.seqname = seqname
        self.source = source
        self.feature = feature
        self.start = long(start)
        self.end = long(end)
        try:
            self.score = float(score)
        except ValueError:
            self.score = None
        self.strand = strand
        try:
            self.frame = int(frame)
        except ValueError:
            self.frame = None
        self.attr = attr

    def __str__(self):
        return "\t".join( [
            self.seqname,
            self.source,
            self.feature,
            str(self.start),
            str(self.end),
            '.' if self.score is None else str(self.score) ,
            self.strand,
            '.' if self.frame is None else str(self.frame),
            str(self.attr)
            ])

    def __repr__(self):
        return "%s:%s:%s-%s" % (self.feature, self.seqname, self.start, self.end)

    def gene_id(self):
        return self.attr['gene_id']

    def get_type(self):
        return self.feature

    def __getitem__(self, name):
        return self.attr[name]



def gtfParse(handle, source_filter=None, feature_filter=None):
    for line in handle:
        if not line.startswith("#"):
            row = line.split("\t")
            attr_str = row[8]
            attr = {}
            for s in attr_str.split("; "):
                res = reAttr.search(s)
                if res:
                    attr[res.group(1)] = res.group(2)
            if source_filter is None or source_filter==row[1]:
                if feature_filter is None or feature_filter==row[2]:
                    yield GTFLine(row[0], row[1], row[2], row[3], row[4], row[5], row[6], row[7], attr)

class GTFGene:
    def __init__(self, gene_id):
        self.elements = {}
        self.gene_id = gene_id
        self.start = None
        self.end = None
        self.seqname = None

    def append(self, elem):
        if self.start is None or elem.start < self.start:
            self.start = elem.start
        if self.end is None or elem.end > self.end:
            self.end = elem.end
        self.seqname = elem.seqname
        t = elem.get_type()
        if t not in self.elements:
            self.elements[t] = []
        self.elements[t].append(elem)

    def gene_total_length(self):
        if 'start_codon' in self.elements and 'stop_codon' in self.elements:
            end = self.elements['stop_codon'][0].end
            start = self.elements['start_codon'][0].start
            return max(start,end) - min(start,end)

    def exon_total_length(self):
        total = 0
        if 'exon' in self.elements:
            for ex in self.elements['exon']:
                total += ex.end - ex.start
        return total



class GTFMap:
    def __init__(self):
        self.gene_map = {}
        self.chrom_map = {}

    def read(self, handle, source_filter=None, feature_filter=None):
        for line in gtfParse(handle, source_filter=source_filter, feature_filter=feature_filter):
            g = line.attr['gene_name'] #.gene_id()

            if line.seqname not in self.chrom_map:
                self.chrom_map[line.seqname] = []

            if g not in self.gene_map:
                gene = GTFGene(g)
                self.gene_map[g] = gene
                self.chrom_map[line.seqname].append(gene)

            self.gene_map[g].append(line)

    def __iter__(self):
        for k in self.gene_map:
            yield k

    def __getitem__(self, n):
        return self.gene_map[n]

    def keys(self):
        return self.gene_map.keys()

    def has_item(self, gene, item):
        if gene not in self.gene_map:
            return False
        if item not in self.gene_map[gene]:
            return False
        return True

    def has_chrom(self, chrome):
        return chrome in self.chrom_map

    def get_chrom(self, chrome):
        return self.chrom_map[chrome]

    def get_item_len(self, gene, item):
        total = 0
        for i in self.gene_map[gene][item]:
            total += i.end - i.start
        return total

    def get_item_span(self, gene, item):
        start = None
        end = None
        for i in self.gene_map[gene][item]:
            if start is None:
                start = i.start
            else:
                start = min(start, i.start)
            if end is None:
                end = i.end
            else:
                end = max(end, i.end)
        return start, end

SEG_SAMPLE   = 0
SEG_SEQ      = 1
SEG_START    = 2
SEG_STOP     = 3
SEG_PROBES   = 4
SEG_VALUE    = 5

def gid_biosample(sample):
    return 'biosample:CCLE:' + sample

def gid_cnacallset_set(sample):
    return "cnaCallSet:biosample:CCLE:" + sample

def parse(seg_file, gtf_file, out_base):
    with gzip.GzipFile(gtf_file) as handle:
        g = GTFMap()
        g.read(handle, source_filter="protein_coding", feature_filter="gene")

    out_handles = {}
    def emit_json(message):
        if message.DESCRIPTOR.full_name not in out_handles:
            out_handles[message.DESCRIPTOR.full_name] = open("ccle." + message.DESCRIPTOR.full_name + ".json", "w")
        msg = json.loads(json_format.MessageToJson(message))
        out_handles[message.DESCRIPTOR.full_name].write(json.dumps(msg))
        out_handles[message.DESCRIPTOR.full_name].write("\n")

    i_map = {}
    for gene_name in g:
        gene = g[gene_name]
        #if args.fix_chrome and gene.seqname.startswith("chr"):
        #    seqname = gene.seqname[3:]
        #else:
        seqname = gene.seqname
        if seqname not in i_map:
            i_map[seqname] = Intersecter()
        i_map[seqname].add_interval( Interval(gene.start, gene.end, value=gene) )

    callset_set = set()
    with open(seg_file) as handle:
        handle.readline()
        for line in handle:
            segment = cna.CNASegment()
            row = line.split("\t")
            cnacallset_id = gid_cnacallset_set(row[SEG_SAMPLE])
            if cnacallset_id not in callset_set:
                callset = cna.CNACallSet()
                callset.id = cnacallset_id
                callset.bio_sample_id = gid_biosample(row[SEG_SAMPLE])
                emit_json(callset)
                callset_set.add(cnacallset_id)

            chrom = row[SEG_SEQ]
            start = long(float(row[SEG_START]))
            stop = long(float(row[SEG_STOP]))
            segment.reference_name = chrom
            segment.start = start
            segment.end = stop
            segment.value = float(row[SEG_VALUE])
            segment.call_set_id = cnacallset_id
            if chrom in i_map:
                for hit in i_map[chrom].find(start, stop+1):
                    segment.genes.append('gene:' + hit.value.gene_id)
            emit_json(segment)

def message_to_json(message):
    msg = json_format.MessageToDict(message)
    msg['#label'] = message.DESCRIPTOR.name
    return json.dumps(msg)

if __name__ == "__main__":
    seg_file = sys.argv[1]
    gtf_file = sys.argv[2]
    out_base = sys.argv[3]
    parse(seg_file, gtf_file, out_base)
