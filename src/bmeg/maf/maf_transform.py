""" transform a maf file into vertex[variant, allele]   """

import logging
import csv
import gzip

from bmeg.vertex import Biosample
from bmeg.edge import CallsetFor, AlleleIn

import bmeg.maf.allele_harvester as allele_harvester

from concurrent.futures import ThreadPoolExecutor
from concurrent.futures import as_completed
from more_itertools import chunked
from itertools import islice

# center = 2
# ncbi_build = 3
CHROMOSOME = "Chromosome"  # 4
START = ["Start_Position", "Start_position"]  # 5
END = ["End_Position", "End_position"]  # 6
# strand = 7
VARIANT_TYPE = "Variant_Type"  # 9
REFERENCE_ALLELE = "Reference_Allele"  # 10
TUMOR_ALLELE = "Tumor_Seq_Allele1"  # 11
tumor_allele2 = "Tumor_Seq_Allele2"  # 12
# annotation_transcript = "Annotation_Transcript"  # 14

TUMOR_SAMPLE_BARCODE = "Tumor_Sample_Barcode"  # 15
NORMAL_SAMPLE_BARCODE = "Matched_Norm_Sample_Barcode"  # 16

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

FEATURE = "Feature"  # 48
FEATURE_TYPE = "Feature_type"  # 49
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


class MAFTransformer():
    def read_maf(self, mafpath, gz, skip=0, harvest=True):
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

    def allele_call_maker(self, allele, line=None):
        """ override, create call from line """
        pass

    def barcode_to_sampleid(self, barcode):  # pragma nocover
        """ override, decode barcode """
        return barcode

    def callset_maker(self, allele, source, centerCol, method, line):  # noqa pragma nocover
        """ override, create callset from line """
        print('override me')
        pass

    def create_gene_gid(self, line):  # pragma nocover
        """ override, create gene_gid from line """
        pass

    def create_allele_dict(self, line, genome='GRCh37'):
        ''' return properly named allele dictionary, populated from line'''
        # collect CURIES that apply to allele
        annotations = []
        annotations.append('{}:{}'.format(VARIANT_TYPE, line.get(VARIANT_TYPE, None)))
        annotations.append('{}:{}'.format(FEATURE_TYPE, line.get(FEATURE_TYPE, None)))
        annotations.append('{}:{}'.format(FEATURE, line.get(FEATURE, None)))
        annotations.append('{}:{}'.format(dbSNP_RS, line.get(dbSNP_RS, None)))
        return {
            'genome': genome,
            'chromosome': line[CHROMOSOME],
            'start': int(get_value(line, START, None)),
            'end': int(get_value(line, END, None)),
            'reference_bases': line[REFERENCE_ALLELE],
            'alternate_bases': line[TUMOR_ALLELE],
            # ,myvariantinfo: dict
            'annotations': annotations,
        }

    def allele_maker(self, line, harvest, filter):
        """ worker task to create and/or harvest allele from line """
        allele_dict = self.create_allele_dict(line)
        return allele_harvester.harvest(**allele_dict,
                                        harvest=harvest,
                                        filter=filter), line

    def multithreading(self, func, lines, max_workers, harvest, filter):
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

    def maf_convert(self, emitter, mafpath, workers, source='tcga',
                    genome='GRCh37',
                    method='variant', gz=False, centerCol='Center', skip=0,
                    harvest=True, filter=[]):
        """
        emitter -  a way to write output
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
        for allele, line in self.multithreading(self.allele_maker,
                                                self.read_maf(
                                                    mafpath, gz, skip),
                                                max_workers=workers,
                                                harvest=harvest,
                                                filter=filter):
            # if allele was filtered out
            if not allele:
                continue
            # save the allele that was created
            emitter.emit_vertex(allele)
            # create edge between the allele and the callset
            call_tuples, callsets = self.callset_maker(allele, source,
                                                       centerCol,
                                                       method, line)
            # save the calls
            for call_tuple in call_tuples:
                call = call_tuple[0]
                callset_gid = call_tuple[1]
                emitter.emit_edge(call, allele.gid(), callset_gid)
            # many callsets can be created, emit only uniques
            for callset in callsets:
                if callset.gid not in my_callsets_ids:
                    my_callsets_ids.add(callset.gid)
                    emitter.emit_vertex(callset)
                    if callset.normal_biosample_id:
                        emitter.emit_edge(CallsetFor(),
                                          callset.gid(),
                                          Biosample.make_gid(callset.normal_biosample_id)
                                          )
                    emitter.emit_edge(CallsetFor(),
                                      callset.gid(),
                                      Biosample.make_gid(callset.tumor_biosample_id)
                                      )

            # create edge to gene
            gene_gid = self.create_gene_gid(line)
            if gene_gid:
                emitter.emit_edge(AlleleIn(), allele.gid(), gene_gid)
            # log progress
            c += 1
            if c % 1000 == 0:  # pragma nocover
                logging.info('imported {}'.format(c))
        logging.info('imported {}'.format(c))
