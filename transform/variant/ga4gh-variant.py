#! /usr/bin/env python
'''
Authors: Malisa Smith smimal@ohsu.edu, Ryan Spangler spanglry@ohsu.edu

'''

import re
import os
import sys
import csv
import gzip
import json
import string
import logging
import requests
import argparse
from bmeg import variants_pb2
from bmeg import conversions
# from ga4gh import bio_metadata_pb2, variants_pb2, allele_annotations_pb2
from google.protobuf import json_format

########################################

def proto_list_append(message, a):
    v = message.values.add()
    v.string_value = a

url = "https://api.gdc.cancer.gov/cases/"
headers = {'Content-Type': "application/json"}
def barcode_id(barcode, retries=5):
    if retries == 0:
        return barcode
    else:
        try:
            query = {
                "expand": "samples",
                "fields": "samples.sample_id",
                "filters": {
                    'op': 'in',
                    'content': {
                        'field': 'samples.portions.analytes.aliquots.submitter_id',
                        'value': [barcode]
                    }
                }
            }

            r = requests.post(url, headers=headers, data=json.dumps(query))
            res = r.json()
            return res['data']['hits'][0].get('samples', [])[0]['sample_id']
        except:
            return barcode_id(barcode, retries - 1)

class Converter(object):
    def __init__(self, 
        bioPrefix, 
        genomePrefix, 
        callSetPrefix, 
        genePrefix, 
        typeField,
        centerCol,
        geneCol,
        ensembl):
        self.bioPrefix = bioPrefix
        self.genomePrefix = genomePrefix
        self.callSetPrefix = callSetPrefix
        self.typeField = typeField
        self.genePrefix = genePrefix
        self.centerCol = centerCol
        self.geneCol = geneCol
        self.ensembl = ensembl

    def gid_biosample(self, name):
        # if self.bioPrefix is not None:
        #     return '%s:%s' % (self.bioPrefix, name)
        return name

    def gid_call_set(self, name):
        return name
        # return '%s:%s' % (self.callSetPrefix, name)

    def gid_variant(self, genome, chromosome, start, end, reference, alternate):
        return '%s:%s:%s:%s:%s:%s' % (genome, chromosome, start, end, reference, alternate)
        # return '%s:%s:%s:%s:%s:%s' % (self.genomePrefix, chromosome, start, end, reference, alternate)

    def gid_variant_annotation(self, variant_id, annotation_set_id):
        return '%s:%s' % (variant_id, annotation_set_id)
        # return '%s:%s:%s' % (self.genomePrefix, variant_id, annotation_set_id)

    def gid_transcript_effect(self, feature_id, annotation_id, alternate_bases):
        return '%s:%s:%s:%s' % (feature_id, annotation_id, alternate_bases)
        # return '%s:%s:%s:%s' % (self.genomePrefix, feature_id, annotation_id, alternate_bases)

    def gid_gene(self, name):
        return '%s:%s' % (self.genePrefix, name)

def get_value(d, keys, default):
    if isinstance(keys, list):
        for k in keys:
            if k in d:
                return d[k]
    elif keys in d:
        return d[keys]
    return default
        

class MafConverter(Converter):
    #def __init__(self, **kwargs):
    #    super(MafConverter, self).__init__(**kwargs)

    def convert(self, emit, mafpath, source, genome, method, gz):
        center = 2
        ncbi_build = 3
        chromosome = "Chromosome" # 4
        start = ["Start_Position", "Start_position"] # 5
        end = ["End_Position", "End_position"] # 6
        strand = 7
        variant_type = "Variant_Type" # 9
        reference_allele = "Reference_Allele" # 10
        tumor_allele1 = "Tumor_Seq_Allele1" # 11
        tumor_allele2 = "Tumor_Seq_Allele2" # 12
        annotation_transcript = "Annotation_Transcript" # 14
        tumor_sample_barcode = "Tumor_Sample_Barcode" # 15

        normal_sample_barcode = "Matched_Norm_Sample_Barcode" # 16
        normal_allele1 = 17
        normal_allele2 = 18
        verification_status = 23
        validation_status = 24
        mutation_status = 25
        sequencing_phase = 26
        sequence_source = 27
        bam_file = 30

        sequencer = 31
        genome_change = 32

        # Information indices for VariantCallEffect and Gene
        entrez_gene_id = 1
        variant_classification = 8

        logging.info('converting maf: ' + mafpath)

        samples = set()
        variants = set()

        sample_ids = {}

        header = None
        if gz:
            inhandle = gzip.GzipFile(mafpath)
        else:
            inhandle = open(mafpath)
        reader = csv.DictReader(inhandle, delimiter="\t")
        for line in reader:
            alternate_bases = set(line[tumor_allele1])
            if tumor_allele2 in line:
                alternate_bases.add( line[tumor_allele2] )

            variant = variants_pb2.Variant()
            variant.start = long(get_value(line, start, None))
            variant.end = long(get_value(line, end, None))
            variant.reference_name = line[chromosome]
            variant.reference_genome = genome
            variant.reference_bases = line[reference_allele]
            #if self.centerCol in line:
            #    for c in line[self.centerCol].split("|"):
                    #proto_list_append(variant.attributes["center"], c)
            #        variant.attributes.get_or_create_list("center").append(c.replace("*", ""))
            #proto_list_append(variant.info["variant_type"], line[variant_type])

            barcode = line[tumor_sample_barcode]
            sample = barcode
            if source == 'tcga':
                if not barcode in sample_ids:
                    sample_ids[barcode] = barcode_id(barcode)
                sample = sample_ids[barcode]

            sample_callsets = []
            if self.centerCol in line:
                for c in line[self.centerCol].split("|"):
                    center = c.replace("*", "")
                    callset_id = "%s:%s" % (sample, center)
                    sample_callsets.append((callset_id, center))
                    call = variant.calls.add()
                    # call.call_set_id = self.gid_call_set(callset_id)
                    call.method = center
                    call.source = source
                    call.biosample_id = sample
            else:
                sample_callsets = [ (sample, None) ]
                call = variant.calls.add()
                # call.call_set_id = self.gid_call_set(sample)
                call.method = method
                call.source = source
                call.biosample_id = sample

            for a in alternate_bases:
                clone = variants_pb2.Variant()
                clone.CopyFrom(variant)
                clone.alternate_bases.append(a)
                variant_id = self.gid_variant(
                    genome,
                    line[chromosome],
                    long(get_value(line, start, None)),
                    long(get_value(line, end, None)),
                    line[reference_allele],
                    a
                )
                clone.id = variant_id
                emit(clone)

                annotation_id = variant_id 
                annotation = variants_pb2.VariantAnnotation()
                annotation.id = annotation_id
                annotation.variant_id = variant_id
                
                feature_id = line[self.geneCol] # self.gid_gene(line[self.geneCol])
                ensembl_id = conversions.hugo_ensembl(feature_id)
                feature_id = feature_id if len(ensembl_id) == 0 else ensembl_id
                
                transcript_effect = annotation.transcript_effects.add()
                transcript_effect.alternate_bases = a
                transcript_effect.feature_id = feature_id
                # transcript_effect.id = self.gid_transcript_effect(
                #     feature_id.replace(self.genePrefix + ":", ""),
                #     annotation_id.replace(self.genomePrefix + ":", ""),
                #     ','.join(alternate_bases))
                
                ontology = transcript_effect.effects.add()
                effect = line[variant_type]
                ontology.term = effect
                #annotations.append(annotation)
                emit(annotation)


            for sc, center in sample_callsets:
                if sc not in samples:
                    samples.add(sc)
                    callset = variants_pb2.CallSet()
                    callset.name = barcode
                    callset.source = source
                    callset.biosample_id = sample # self.gid_biosample(sample)
                    callset.method = center or method
                    # callset.id = self.gid_call_set(sc)
                    # callset.attributes["center"] = center

                    emit(callset)
                    
            # for v in variants:
            #     emit(v)


        inhandle.close()

def convert_to_profobuf(maf, vcf, gz, out, multi, format, ensembl, source, genome, method, bioPrefix, genomePrefix, callSetPrefix, genePrefix, typeField, centerCol, geneCol):
    out_handles = {}
    if maf:
        m = MafConverter(
            bioPrefix=bioPrefix,
            genomePrefix=genomePrefix,
            callSetPrefix=callSetPrefix,
            genePrefix=genePrefix,
            typeField=typeField,
            centerCol=centerCol,
            geneCol=geneCol,
            ensembl=ensembl,
        )
        def emit_json_single(message):
            if 'main' not in out_handles:
                out_handles['main'] = open(out, 'w')
            msg = json.loads(json_format.MessageToJson(message))
            msg[typeField] = message.DESCRIPTOR.full_name
            out_handles['main'].write(json.dumps(msg))
            out_handles['main'].write("\n")
        def emit_json_multi(message):
            if message.DESCRIPTOR.full_name not in out_handles:
                out_handles[message.DESCRIPTOR.full_name] = open(multi + '.' + message.DESCRIPTOR.full_name + '.json', 'w')
            msg = json.loads(json_format.MessageToJson(message))
            out_handles[message.DESCRIPTOR.full_name].write(json.dumps(msg))
            out_handles[message.DESCRIPTOR.full_name].write('\n')

        emit = emit_json_single if out else emit_json_multi
        if os.path.isdir(maf):
            for path in os.listdir(maf):
                if path[-4:] == '.maf':
                    m.convert(emit, maf + '/' + path, source, genome, method)
        else:
            m.convert(emit, maf, source, genome, method, gz)

    for handle in out_handles.values():
        handle.close()

def parse_args(args):
    # We don't need the first argument, which is the program name
    args = args[1:]

    # Construct the parser
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    # Now add all the options to it
    parser.add_argument('--gz', action="store_true", default=False, help='Path to the maf you want to import')
    parser.add_argument('--maf', type=str, help='Path to the maf you want to import')
    parser.add_argument('--vcf', type=str, help='Path to the vcf you want to import')
    parser.add_argument('--ensembl', type=str, help='Ensembl Map (see ftp://ftp.ensembl.org/pub/release-90/tsv/homo_sapiens/Homo_sapiens.GRCh38.90.entrez.tsv.gz )')
    parser.add_argument('--out', type=str, help='Path to output file (.json or .pbf_ext)')
    parser.add_argument('--multi', type=str, help='Path to output file prefix (.json or .pbf_ext)')
    parser.add_argument('--source', default='mc3', help='source of the sample')
    parser.add_argument('--genome', default='GRCh37', help='source of the sample')
    parser.add_argument('--method', default='variant', help='method used to process the variants')
    parser.add_argument('--bioPrefix', default=None)
    parser.add_argument('--genomePrefix', default='grch37')
    parser.add_argument('--callSetPrefix', default='callSet')
    parser.add_argument('--typeField', default='#label')
    parser.add_argument('--genePrefix', default='gene')
    parser.add_argument('--center', dest="centerCol", default='Center', help='caller field')
    parser.add_argument('--gene', dest="geneCol", default="Entrez_Gene_Id")
    parser.add_argument('--format', type=str, default='json', help='Format of output: json or pbf (binary)')
    return parser.parse_args(args)

if __name__ == '__main__':
    options = parse_args(sys.argv)
    convert_to_profobuf(**vars(options))
