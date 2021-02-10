"""
dedupe and harmonize alleles
"""
import logging
import glob
import sys
import subprocess
import os.path
import pandas

from bmeg.util.cli import default_argument_parser
from bmeg.util.logging import default_logging
from bmeg.emitter import new_emitter
from bmeg import (Gene, Transcript, Protein, PfamFamily,
                  Allele_Gene_Gene, Allele_Transcript_Transcript,
                  Allele_Protein_Protein, Allele_PfamFamily_PfamFamily)

from bmeg.maf import make_allele


def valid_filenames(path):
    """ filter out any we do not want to process """
    filenames = []
    print("Search: %s" % (path))
    files = glob.glob(path, recursive=True)
    for filename in files:
        # myvariant download **DEPRECIATED**
        if 'myvariant' in filename:
            continue
        # pseudo alleles **DEPRECIATED**
        elif 'MinimalAllele' in filename:
            continue
        # our own output
        elif 'normalized.' in filename:
            continue
        else:
            filenames.append(filename)
    return filenames


def create_minimal_maf(path, minimal_maf):
    """ dedup alleles and create a minimal MAF file """
    if os.path.isfile(minimal_maf):
        logging.info('FOUND: {}'.format(minimal_maf))
        return
    try:
        files = valid_filenames(path)
        assert len(files) > 0
        files = ' '.join(files)
        logging.info('FILES: {}'.format(files))
        cmd = "echo 'Chromosome\tStart_Position\tReference_Allele\tTumor_Seq_Allele2\tTumor_Sample_Barcode' > {}".format(minimal_maf)
        logging.info('RUNNING: {}'.format(cmd))
        subprocess.check_output(cmd, shell=True)
        cmd = "zcat %s | jq -cr '[.data.chromosome, .data.start, .data.reference_bases, .data.alternate_bases, \"STUBID\"] | @tsv' | grep -v -E '^M|^X|^Y' | sort -n -k 1,1 -k 2,2 | uniq  >> %s" % (files, minimal_maf)
        logging.info('RUNNING: {}'.format(cmd))
        subprocess.check_output(cmd, shell=True)
        cmd = "zcat %s | jq -cr '[.data.chromosome, .data.start, .data.reference_bases, .data.alternate_bases, \"STUBID\"] | @tsv' | grep -E '^M|^X|^Y' | sort -k 1,1 -k2,2 | sed 's/^M/MT/g' | uniq >> %s" % (files, minimal_maf)
        logging.info('RUNNING: {}'.format(cmd))
        subprocess.check_output(cmd, shell=True)
        return
    except subprocess.CalledProcessError as sort_error:
        raise ValueError("ERROR [{}]: {}".format(sort_error.returncode, sort_error.output))


def run_maf2maf(minimal_maf, annotated_maf):
    if os.path.isfile(annotated_maf):
        logging.info('FOUND: {}'.format(annotated_maf))
        return
    try:
        cmd = "bash transform/allele/run_maf2maf.sh {} {}".format(minimal_maf, annotated_maf)
        logging.info('RUNNING: {}'.format(cmd))
        subprocess.check_output(cmd, shell=True)
        return
    except subprocess.CalledProcessError as maf2maf_error:
        raise ValueError("ERROR [{}]: {}".format(maf2maf_error.returncode, maf2maf_error.output))


def transform(output_dir='pre-outputs',
              vertex_filename_pattern='*/*Allele.Vertex.json.gz',
              minimal_maf_file='pre-outputs/allele-maf/minimal_alleles.maf',
              annotated_maf_file='pre-outputs/allele-maf/annotated_alleles.maf',
              emitter_directory='allele',
              emitter_prefix='normalized',
              emitter_name='json'):
    """
    dedup alleles & harmonize annotations.
    """

    emitter = new_emitter(name=emitter_name, directory=emitter_directory, prefix=emitter_prefix)

    # Step one:
    # sort and dedupe the outputs/**/*.Allele.Vertex.json.gz files
    if not os.path.isfile(minimal_maf_file):
        path = '{}/{}'.format(output_dir, vertex_filename_pattern)
        create_minimal_maf(path, minimal_maf_file)
    else:
        logging.info('using existing minimal maf file: %s' % minimal_maf_file)

    # Step two:
    # run maf2maf on the minimal MAF file from step one
    if not os.path.isfile(annotated_maf_file):
        run_maf2maf(minimal_maf_file, annotated_maf_file)
    else:
        logging.info('using existing annotated maf file: %s' % annotated_maf_file)

    # Step three:
    # emit all alleles and corresponding edges
    logging.info('emitting annotated alleles')
    ac = ec = 0
    maf = pandas.read_csv(annotated_maf_file, sep='\t', comment='#', dtype=str)
    for index, line in maf.iterrows():
        line = line.dropna().to_dict()

        try:
            allele = make_allele(line)
            emitter.emit_vertex(allele)
            ac += 1
        except Exception as e:
            logging.error("emit allele: %s" % str(e))

        if allele.ensembl_gene:
            try:
                emitter.emit_edge(
                    Allele_Gene_Gene(
                        from_gid=allele.gid(),
                        to_gid=Gene.make_gid(allele.ensembl_gene)
                    ),
                    emit_backref=True
                )
                ec += 1
            except Exception as e:
                logging.error("emit allele -> gene edge: %s" % str(e))

        if allele.ensembl_transcript:
            try:
                emitter.emit_edge(
                    Allele_Transcript_Transcript(
                        from_gid=allele.gid(),
                        to_gid=Transcript.make_gid(allele.ensembl_transcript)
                    ),
                    emit_backref=True
                )
                ec += 1
            except Exception as e:
                logging.error("emit allele -> transcript edge: %s" % str(e))

        if allele.ensembl_protein:
            try:
                emitter.emit_edge(
                    Allele_Protein_Protein(
                        from_gid=allele.gid(),
                        to_gid=Protein.make_gid(allele.ensembl_protein)
                    ),
                    emit_backref=True
                )
                ec += 1
            except Exception as e:
                logging.error("emit allele -> protein edge: %s" % str(e))

        if 'Pfam_domain' in line.get('DOMAINS', ''):
            try:
                domains = [x.split(':')[1] for x in line.get('DOMAINS', '').split(',') if 'Pfam_domain' in x]
                for d in domains:
                    emitter.emit_edge(
                        Allele_PfamFamily_PfamFamily(
                            from_gid=allele.gid(),
                            to_gid=PfamFamily.make_gid(d)
                        ),
                        emit_backref=True
                    )
                    ec += 1
            except Exception as e:
                logging.error("emit allele -> PFAM edge: %s" % str(e))

        if ac % 10000 == 0:
            logging.info('alleles emitted: {}'.format(ac))
            logging.info('edges emitted: {}'.format(ec))

    logging.info('total alleles emitted: {}'.format(ac))
    logging.info('total edges emitted: {}'.format(ec))

    emitter.close()


def main():  # pragma: no cover
    parser = default_argument_parser(emitter_directory_default='allele',
                                     emitter_prefix_default='normalized')
    parser.add_argument('--allele-output-dir', type=str,
                        help='Path to the directory containing */*.Allele.Vertex.json.gz',
                        default='pre-outputs')
    parser.add_argument('--vertex-file-pattern', type=str,
                        help='vertex file pattern to glob for',
                        default='*/*Allele.Vertex.json.gz')
    parser.add_argument('--minimal-maf', type=str,
                        help='intermediate minimal maf file',
                        default='pre-outputs/allele-maf/minimal_alleles.maf')
    parser.add_argument('--annotated-maf', type=str,
                        help='intermediate annotated maf file',
                        default='pre-outputs/allele-maf/annotated_alleles.maf')
    parser.add_argument('--leave-intermediate-files', dest='delete_intermediate', action='store_false',
                        help='do not delete minimal or annotated mafs prior to running')
    parser.set_defaults(delete_intermediate=True)

    # We don't need the first argument, which is the program name
    options = parser.parse_args(sys.argv[1:])
    default_logging(options.loglevel)

    if options.delete_intermediate and os.path.isfile(options.minimal_maf):
        logging.info('deleting {}'.format(options.minimal_maf))
        os.remove(options.minimal_maf)

    if options.delete_intermediate and os.path.isfile(options.annotated_maf):
        logging.info('deleting {}'.format(options.annotated_maf))
        os.remove(options.annotated_maf)

    transform(output_dir=options.allele_output_dir,
              vertex_filename_pattern=options.vertex_file_pattern,
              minimal_maf_file=options.minimal_maf,
              annotated_maf_file=options.annotated_maf,
              emitter_directory=options.emitter_directory,
              emitter_prefix=options.emitter_prefix,
              emitter_name=options.emitter)


if __name__ == '__main__':  # pragma: no cover
    main()
