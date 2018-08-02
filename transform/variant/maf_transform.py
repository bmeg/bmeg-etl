#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" transform a maf file into vertexs[variant, allele]   """

import logging
import csv
import gzip
import sys

from bmeg.models.vertex_models import Biosample, Callset, Gene
from bmeg.models.edge_models import AlleleCall, CallsetFor, AlleleIn
from bmeg.models.emitter import Emitter
from bmeg.util.cli import default_argument_parser
from bmeg.util.logging import default_logging

import allele_harvester

from concurrent.futures import ThreadPoolExecutor
from concurrent.futures import as_completed
from more_itertools import chunked
from itertools import islice

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


def get_value(d, keys, default):
    """ utility get value from list"""
    if isinstance(keys, list):
        for k in keys:
            if k in d:
                return d[k]
    elif keys in d:
        return d[keys]
    return default


def _read_maf(mafpath, gz, skip=0, harvest=True):
    """ generator for each line in maf """
    if gz or 'gz' in mafpath:
        inhandle = gzip.open(mafpath, mode='rt')
    else:
        inhandle = open(mafpath)
    reader = csv.DictReader(inhandle, delimiter="\t")
    logging.info('skipping: {}'.format(skip))
    reader = islice(reader, skip, None)
    for line in reader:
        yield line
    inhandle.close()


def _allele_call_maker(allele, callset, line=None):
    """ create call from line """
    keys = ['t_depth', 't_ref_count', 't_alt_count', 'n_depth',
            'n_ref_count', 'n_alt_count', 'FILTER',
            'Match_Norm_Seq_Allele1',	'Match_Norm_Seq_Allele2',
            'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2',
            ]
    info = {}
    for k in keys:
        info[k] = get_value(line, k, None)
    return AlleleCall(allele.gid, callset.gid, info)


def _tcga_aliquot2sample_barcode(barcode):
    """ create tcga sample barcode """
    return "-".join(barcode.split("-")[0:4])


def _callset_maker(allele, source, centerCol, method, line):
    """ create callset from line """
    barcode = line[tumor_sample_barcode]
    sample = barcode
    if source == 'tcga':
        sample = _tcga_aliquot2sample_barcode(barcode)
    sample = Biosample.make_gid(sample)
    sample_callsets = []
    sample_calls = []
    if centerCol in line:
        for c in line[centerCol].split("|"):
            center = c.replace("*", "")
            # callset_id = "%s:%s" % (sample, center)
            callset = Callset(sample,
                              Biosample.make_gid(
                                line[normal_sample_barcode]), center)
            sample_callsets.append(callset)
            sample_calls.append(
                _allele_call_maker(allele, callset, line))
    else:
        callset = Callset(sample,
                          Biosample.make_gid(
                            line[normal_sample_barcode]), method)
        sample_callsets.append(callset)
        sample_calls.append(
            _allele_call_maker(allele, callset, line))
    return sample_calls, sample_callsets


def _allele_dict(line, genome='GRCh37'):
    ''' return properly named allele dictionary, populated form line'''
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

    return {
        'genome': genome,
        'chromosome': line[chromosome],
        'start': int(get_value(line, start, None)),
        'end': int(get_value(line, end, None)),
        'reference_bases': line[reference_allele],
        'alternate_bases': line[tumor_allele2],
        # ,myvariantinfo: dict
        'annotations': annotations,
    }


def _allele_maker(line, harvest, filter):
    """ worker task to create and/or harvest allele from line """
    allele_dict = _allele_dict(line)
    return allele_harvester.harvest(**allele_dict,
                                    harvest=harvest,
                                    filter=filter), line
    # if harvest:
    #     return allele_harvester.harvest(**allele_dict), line
    # return allele_harvester.create(**allele_dict), line


def _multithreading(func, lines, max_workers, harvest, filter):
    """
    Create a thread pool and create alleles
    """
    # limit queue size == number of workers
    # see https://stackoverflow.com/questions/48263704/threadpoolexecutor-how-to-limit-the-queue-maxsize  # noqa

    for chunk in chunked(lines, max_workers):
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = (executor.submit(func, line, harvest, filter) for line in chunk)  # noqa
            for future in as_completed(futures):
                yield future.result()


def _allele_to_gene(allele, line):
    """
    create edge to gene
    """
    ensembl_id = line.get('Gene', None)
    symbol = line.get('Hugo_Symbol', None)
    gene_gid = Gene.make_gid(ensembl_id=ensembl_id)
    return AlleleIn(allele.gid, gene_gid)


def maf_convert(emit, mafpath, workers, source='tcga', genome='GRCh37',
                method='variant', gz=False, centerCol='Center', skip=0,
                harvest=True, filter=[]):
    """
    emit -  a way to write output
    mafpath - a file to read
    source - special handling if 'tcga'
    genome - reference_genome e.g. GRCh37
    method - call's method
    gz - is mafpath a gz
    skip - ignore first N lines
    """

    logging.info('converting maf: ' + mafpath)
    my_callsets_ids = set()
    c = skip
    for allele, line in _multithreading(_allele_maker,
                                        _read_maf(mafpath, gz, skip),
                                        max_workers=workers,
                                        harvest=harvest,
                                        filter=filter):
        # if allele was filtered out
        if not allele:
            continue
        # save the allele that was created
        emit(allele)
        # create edge between the allele and the callset
        calls, callsets = _callset_maker(allele, source, centerCol,
                                         method, line)
        # save the calls
        for call in calls:
            emit(call)
        # many callsets can be created, emit only uniques
        for callset in callsets:
            if callset.gid not in my_callsets_ids:
                my_callsets_ids.add(callset.gid)
                emit(callset)
                emit(CallsetFor(callset.gid, callset.normal_biosample_id))
                emit(CallsetFor(callset.gid, callset.tumor_biosample_id))

        # create edge to gene
        emit(_allele_to_gene(allele, line))
        # log progress
        c += 1
        if c % 1000 == 0:  # pragma nocover
            logging.info('imported {}'.format(c))
    logging.info('imported {}'.format(c))


def convert(mafpath, prefix, workers=5, skip=0, harvest=True, filter=[]):
    """ entry point """
    emitter = Emitter(prefix=prefix)
    maf_convert(emit=emitter.emit, mafpath=mafpath, workers=workers, skip=skip,
                harvest=harvest, filter=filter)


if __name__ == '__main__':  # pragma: no cover
    parser = default_argument_parser()
    parser.add_argument('--maf_file', type=str,
                        help='Path to the maf you want to import')
    parser.add_argument(
        '--workers', type=int,
        help="multithread harvest from myvariant.info",
        default=5
    )
    parser.add_argument(
        '--skip', type=int,
        help="skip first N lines in MAF",
        default=0
    )
    parser.add_argument('--filter', type=str,
                        help='Path of already harvested Allele gids')
    harvest = parser.add_mutually_exclusive_group(required=False)
    harvest.add_argument('--harvest', dest='harvest',
                         help="get myvariantinfo", action='store_true')
    harvest.add_argument('--no-harvest', dest='harvest',
                         help="do not get myvariantinfo", action='store_false')
    parser.set_defaults(harvest=True)

    #
    # parser.add_argument('--harvest', dest='harvest', action='store_true',
    #                     help="retrieve external data",
    #                     default=True)

    # We don't need the first argument, which is the program name
    options = parser.parse_args(sys.argv[1:])
    default_logging(options.loglevel)

    # ids to skip
    filter = {}
    if options.filter:
        with open(options.filter) as f:
            for line in f:
                line = line.strip()
                filter[line] = None

    convert(mafpath=options.maf_file, prefix=options.prefix,
            workers=options.workers, skip=options.skip,
            harvest=options.harvest, filter=filter)
