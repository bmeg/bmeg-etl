#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" transform a maf file into vertexs[variant, allele]   """

import argparse
import logging
import csv
import gzip
import sys

from bmeg.models.vertex_models import Biosample, Callset
from bmeg.models.edge_models import AlleleCall, CallsetFor
from bmeg.models.emitter import Emitter

import allele_harvester

from concurrent.futures import ThreadPoolExecutor
from concurrent.futures import as_completed
from itertools import islice


def get_value(d, keys, default):
    """ utility get value from list"""
    if isinstance(keys, list):
        for k in keys:
            if k in d:
                return d[k]
    elif keys in d:
        return d[keys]
    return default


def _multithreading(func, lines, max_workers=5):
    """
    Create a thread pool and create alleles
    """
    # limit queue size == number of workers
    # see https://stackoverflow.com/questions/48263704/threadpoolexecutor-how-to-limit-the-queue-maxsize  # noqa

    def chunked_iterable(iterable, chunk_size):
        it = iter(iterable)
        while True:
            chunk = tuple(islice(it, chunk_size))
            if not chunk:
                break
            yield chunk

    for chunk in chunked_iterable(lines, max_workers):
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = (executor.submit(func, line) for line in chunk)
            for future in as_completed(futures):
                yield future.result()


def _read_maf(mafpath, gz):
    """ generator for each line in maf """
    if gz or 'gz' in mafpath:
        inhandle = gzip.open(mafpath, mode='rt')
    else:
        inhandle = open(mafpath)
    reader = csv.DictReader(inhandle, delimiter="\t")
    for line in reader:
        yield line
    inhandle.close()


class MafConverter():

    def convert(self, emit, mafpath, source, genome='GRCh37',
                method='variant', gz=False, centerCol='Center'):
        """
        emit -  a way to write output
        mafpath - a file to read
        source - special handling if 'tcga'
        genome - reference_genome e.g. GRCh37
        method - call's method
        gz - is mafpath a gz
        """
        # center = 2
        # ncbi_build = 3
        chromosome = "Chromosome"  # 4
        start = ["Start_Position", "Start_position"]  # 5
        end = ["End_Position", "End_position"]  # 6
        # strand = 7
        variant_type = "Variant_Type"  # 9
        reference_allele = "Reference_Allele"  # 10
        # tumor_allele1 = "Tumor_Seq_Allele1"  # 11
        tumor_allele2 = "Tumor_Seq_Allele2"  # 12
        # annotation_transcript = "Annotation_Transcript"  # 14

        tumor_sample_barcode = "Tumor_Sample_Barcode"  # 15
        normal_sample_barcode = "Matched_Norm_Sample_Barcode"  # 16

        # normal_allele1 = 17
        # normal_allele2 = 18
        # verification_status = 23
        # validation_status = 24
        # mutation_status = 25
        # sequencing_phase = 26
        # sequence_source = 27
        # bam_file = 30

        # sequencer = 31
        # genome_change = 32

        # # Information indices for VariantCallEffect and Gene
        # entrez_gene_id = 1
        # variant_classification = 8

        feature = "Feature"  # 48
        feature_type = "Feature_type"  # 49
        dbSNP_RS = "dbSNP_RS"  # 13

        logging.info('converting maf: ' + mafpath)

        sample_ids = {}

        def _allele_call_maker(allele, callset, line=None):
            """ create call from line """
            keys = ['t_depth', 't_ref_count', 't_alt_count', 'n_depth',
                    'n_ref_count', 'n_alt_count']
            info = {}
            for k in keys:
                info[k] = get_value(line, k, None)
            return AlleleCall(allele.gid, callset.gid, info)

        def _callset_maker(allele, line):
            """ create callset from line """
            barcode = line[tumor_sample_barcode]
            sample = barcode
            if source == 'tcga':
                if barcode not in sample_ids:
                    sample_ids[barcode] = _tcga_aliquot2sample_barcode(barcode)
                sample = sample_ids[barcode]
            sample = Biosample.make_gid(sample)
            _sample_callsets = []
            _sample_calls = []
            if centerCol in line:
                for c in line[centerCol].split("|"):
                    center = c.replace("*", "")
                    # callset_id = "%s:%s" % (sample, center)
                    callset = Callset(sample,
                                      Biosample.make_gid(
                                        line[normal_sample_barcode]), center)
                    _sample_callsets.append(callset)
                    _sample_calls.append(
                        _allele_call_maker(allele, callset, line))
            else:
                callset = Callset(sample,
                                  Biosample.make_gid(
                                    line[normal_sample_barcode]), method)
                _sample_callsets.append(callset)
                _sample_calls.append(
                    _allele_call_maker(allele, callset, line))
            return _sample_calls, _sample_callsets

        def _allele_maker(line):
            """ worker task to create allele from line """
            # collect CURIES that apply to allele
            annotations = []
            annotations.append('{}:{}'.format(
                variant_type, line[variant_type]))
            annotations.append('{}:{}'.format(
                feature_type, line[feature_type]))
            annotations.append('{}:{}'.format(
                feature, line[feature]))
            annotations.append('{}:{}'.format(
                dbSNP_RS, line[dbSNP_RS]))

            allele_dict = {
                'genome': genome,
                'chromosome': line[chromosome],
                'start': int(get_value(line, start, None)),
                'end': int(get_value(line, end, None)),
                'reference_bases': line[reference_allele],
                'alternate_bases': line[tumor_allele2],
                # ,myvariantinfo: dict
                'annotations': annotations,
            }
            return allele_harvester.harvest(**allele_dict), line

        def _tcga_aliquot2sample_barcode(barcode):
            """ create tcga sample barcode """
            return "-".join(barcode.split("-")[0:4])

        my_callsets_ids = set()
        my_callsets = []

        for allele, line in _multithreading(_allele_maker,
                                            _read_maf(mafpath, gz)):
            # save the allele that was created
            emit(allele)

            # create edge between the allele and the callset
            calls, _callsets = _callset_maker(allele, line)
            # many callsets can be created, save until end of maf processing
            # to emit only uniques
            for callset in _callsets:
                if callset.gid not in my_callsets_ids:
                    my_callsets_ids.add(callset.gid)
                    my_callsets.append(callset)
            # emit calls now
            for call in calls:
                emit(call)

        for callset in my_callsets:
            emit(callset)
            emit(CallsetFor(callset.gid, callset.normal_biosample_id))
            emit(CallsetFor(callset.gid, callset.tumor_biosample_id))


def convert(maf_file, prefix):
    converter = MafConverter()
    emitter = Emitter(prefix=prefix)
    converter.convert(emit=emitter.emit, mafpath=maf_file, source='tcga')


def parse_args(args):  # pragma: no cover
    # We don't need the first argument, which is the program name
    args = args[1:]

    # Construct the parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    # Now add all the options to it
    # parser.add_argument('--gz', action="store_true", default=False,
    #                     help='Path to the maf you want to import')
    parser.add_argument('--maf_file', type=str,
                        help='Path to the maf you want to import')
    parser.add_argument('--prefix', type=str,
                        help='Path prefix for output files')
    return parser.parse_args(args)


if __name__ == '__main__':  # pragma: no cover
    options = parse_args(sys.argv)
    fmt = '%(asctime)s - %(levelname)s - %(threadName)s - %(message)s'
    logging.basicConfig(level=logging.INFO, format=fmt)
    convert(**vars(options))
