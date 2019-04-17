""" transform a maf file into vertex[variant, allele]   """

import logging
import csv
import gzip
import sys

from bmeg import Allele, Aliquot
from bmeg.gid import make_allele_gid
#from bmeg.edge import HasCallset, AlleleIn
from bmeg.emitter import new_emitter
from bmeg.util.cli import default_argument_parser
from bmeg.util.logging import default_logging


from concurrent.futures import ThreadPoolExecutor
from concurrent.futures import as_completed
from more_itertools import chunked
from itertools import islice


STANDARD_MAF_KEYS = {
    'Hugo_Symbol': 'hugo_symbol',
    'Transcript_ID': 'ensembl_transcript',
    'Variant_Classification': 'effect',
    'Variant_Type': 'type',
    'dbSNP_RS': 'dbSNP_RS'
}

# center = 2
# ncbi_build = 3
CHROMOSOME = "Chromosome"  # 4
START = ["Start_Position", "Start_position"]  # 5
END = ["End_Position", "End_position"]  # 6

STRAND = "Strand"  # strand = 7
VARIANT_TYPE = "Variant_Type"  # 9
REFERENCE_ALLELE = "Reference_Allele"  # 10
tumor_allele1 = "Tumor_Seq_Allele1"  # 11
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
    # override the column used for tumor allele
    TUMOR_ALLELE = tumor_allele1
    # options argument
    DEFAULT_MAF_FILE = None

    def read_maf(self, mafpath, gz, skip=0, harvest=True):
        """ generator for each line in maf """
        if gz or 'gz' in mafpath:
            inhandle = gzip.open(mafpath, mode='rt')
        else:
            inhandle = open(mafpath)
        self.current_path = mafpath
        logging.info('skipping {}'.format(skip))
        inh = islice(inhandle, skip, None)
        reader = csv.DictReader(inh, delimiter="\t")

        for line in reader:
            yield line
        inhandle.close()

    def allele_call_maker(self, allele, line=None):
        """ override, create call from line """
        pass

    def barcode_to_aliquot_id(self, barcode):  # pragma nocover
        """ override, decode barcode """
        return barcode

    def callset_maker(self, allele, source, centerCol, method, line):  # noqa pragma nocover
        """ override, create callset from line """
        logging.error('override me')
        pass

    def create_gene_gid(self, line):  # pragma nocover
        """ override, create gene_gid from line """
        pass

    def create_allele_dict(self, line, genome='GRCh37'):
        ''' return properly named allele dictionary, populated from line'''
        # collect CURIES that apply to allele
        record = {
            'genome': genome,
            'chromosome': line[CHROMOSOME],
            'start': int(get_value(line, START, None)),
            'end': int(get_value(line, END, None)),
            'reference_bases': line[REFERENCE_ALLELE],
            'alternate_bases': line[self.TUMOR_ALLELE],
            'strand': line[STRAND]
        }
        for key, data_key in STANDARD_MAF_KEYS.items():
            value = line.get(key, None)
            if value:
                record[data_key] = value
        return record

    def allele_maker(self, line):
        """ worker task to create and/or harvest allele from line """
        allele_dict = self.create_allele_dict(line)
        allele_id = make_allele_gid(
            genome=allele_dict["genome"],
            chromosome=allele_dict["chromosome"],
            start=allele_dict["start"],
            end=allele_dict["end"],
            ref_bases=allele_dict["reference_bases"],
            alt_bases=allele_dict["alternate_bases"])
        return Allele(allele_id=allele_id, **allele_dict)

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

    def maf_convert(self, emitter, mafpath, source,
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
        e = 0
        for line in self.read_maf(mafpath, gz, skip):
            try:
                allele = self.allele_maker(line)
                # save the allele that was created
                emitter.emit_vertex(allele)
                # create edge between the allele and the callset
                call_tuples, callsets = self.callset_maker(allele, source, centerCol, method, line)
                # save the calls
                for call_tuple in call_tuples:
                    call = call_tuple[0]
                    callset_gid = call_tuple[1]
                    emitter.emit_edge(call, from_gid=callset_gid, to_gid=allele.allele_id)
                # many callsets can be created, emit only uniques
                for callset in callsets:
                    if callset.id not in my_callsets_ids:
                        my_callsets_ids.add(callset.id)
                        emitter.emit_vertex(callset)
                        if callset.normal_aliquot_id:
                            emitter.emit_link("hasCallset",
                                              to_gid=callset.id,
                                              from_gid=callset.normal_aliquot_id #Aliquot.make_gid(callset.normal_aliquot_id),
                                              )
                        emitter.emit_link("hasCallset",
                                          to_gid=callset.id,
                                          from_gid=callset.tumor_aliquot_id #Aliquot.make_gid(callset.tumor_aliquot_id),
                                          )

                # create edge to gene
                gene_gid = self.create_gene_gid(line)
                if gene_gid:
                    emitter.emit_link("alleleIn", allele.allele_id, gene_gid)
            except Exception as exc:
                logging.exception(exc)
                e += 1
                #emitter.emit_vertex(Deadletter(target_label='Allele', data=line))

            # log progress
            c += 1
            if c % 1000 == 0:  # pragma nocover
                logging.info('imported {} errors: {}'.format(c, e))
        logging.info('imported {}'.format(c))


def transform(mafpath, prefix, source, emitter_name='json', skip=0, transformer=MAFTransformer()):
    """ entry point """
    emitter = new_emitter(name=emitter_name, directory=prefix)
    transformer.maf_convert(emitter=emitter, mafpath=mafpath, skip=skip, source=source)
    emitter.close()


def maf_default_argument_parser(transformer):
    """ add our default arguments """
    parser = default_argument_parser(transformer.DEFAULT_PREFIX)
    parser.add_argument('--maf_file', type=str,
                        help='Path to the maf you want to import',
                        default=transformer.DEFAULT_MAF_FILE)
    parser.add_argument(
        '--skip', type=int,
        help="skip first N lines in MAF",
        default=0
    )
    return parser


def main(transformer, maf_file=None):  # pragma: no cover
    parser = maf_default_argument_parser(transformer)
    # We don't need the first argument, which is the program name
    options = parser.parse_args(sys.argv[1:])
    default_logging(options.loglevel)

    transform(mafpath=options.maf_file,
              prefix=options.prefix,
              skip=options.skip,
              source=transformer.SOURCE,
              emitter_name=options.emitter,
              transformer=transformer)


if __name__ == '__main__':  # pragma: no cover
    main(transformer=MAFTransformer)
