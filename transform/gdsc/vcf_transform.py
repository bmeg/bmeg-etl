from glob import glob
from types import SimpleNamespace as SN
import logging
import os

import bmeg.ioutils
from bmeg.util.cli import default_argument_parser
from bmeg.util.logging import default_logging
from bmeg import Aliquot, Project, SomaticCallset, SomaticCallset_Alleles_Allele, SomaticCallset_Aliquots_Aliquot
from bmeg.emitter import new_emitter
from bmeg.vcf import read_vcf, make_minimal_allele

PROJECT_NAME = "GDSC1000"

def parse_genotypes(vcf_line, file_type, sample_name):
    """
    Parse format and genotype fields.

    For pindel files:

    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    ##FORMAT=<ID=PP,Number=1,Type=Integer,Description="Pindel calls on the positive strand">
    ##FORMAT=<ID=NP,Number=1,Type=Integer,Description="Pindel calls on the negative strand">
    ##FORMAT=<ID=PB,Number=1,Type=Integer,Description="BWA calls on the positive strand">
    ##FORMAT=<ID=NB,Number=1,Type=Integer,Description="BWA calls on the negative strand">
    ##FORMAT=<ID=PD,Number=1,Type=Integer,Description="BWA mapped reads on the positive strand">
    ##FORMAT=<ID=ND,Number=1,Type=Integer,Description="BWA mapped reads on the negative strand">
    ##FORMAT=<ID=PR,Number=1,Type=Integer,Description="Total mapped reads on the positive strand">
    ##FORMAT=<ID=NR,Number=1,Type=Integer,Description="Total mapped reads on the negative strand">
    ##FORMAT=<ID=PU,Number=1,Type=Integer,Description="Unique calls on the positive strand">
    ##FORMAT=<ID=NU,Number=1,Type=Integer,Description="Unique calls on the negative strand">

    For caveman files:

    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    ##FORMAT=<ID=AA,Number=1,Type=Integer,Description="Reads presenting a A for this position">
    ##FORMAT=<ID=CA,Number=1,Type=Integer,Description="Reads presenting a C for this position">
    ##FORMAT=<ID=GA,Number=1,Type=Integer,Description="Reads presenting a G for this position">
    ##FORMAT=<ID=TA,Number=1,Type=Integer,Description="Reads presenting a T for this position">
    ##FORMAT=<ID=PM,Number=1,Type=Float,Description="Proportion of mut allele">
    """
    gt_format = vcf_line['FORMAT'].split(':')
    gt_info = vcf_line[sample_name].split(':')
    geno = dict(zip(gt_format, gt_info))
    if file_type == 'pindel':
        p = int(geno.get('PP')) if geno.get('PP') else None
        n = int(geno.get('NP')) if geno.get('NP') else None
        d = p + n if p and n else None
        stats = {
            'depth_of_coverage': d,
            'ref_depth': None,
            'alt_depth': None
        }
    elif file_type == 'caveman':
        allele_map = {
            'A': 'AA',
            'C': 'CA',
            'G': 'GA',
            'T': 'TA',
        }
        r = geno.get(allele_map.get(vcf_line.get('REF')))
        r = int(r) if r else None
        a = geno.get(allele_map.get(vcf_line.get('ALT')))
        a = int(a) if a else None
        d = a + r if a and r else None
        stats = {
            'depth_of_coverage': d,
            'ref_depth': r,
            'alt_depth': a
        }
    else:
        raise ValueError("unknown file_type: {}".format(file_type))

    return SN(**stats)


def make_variant_call_data(vcf_line, method):
    tumor_geno = parse_genotypes(vcf_line, method, 'TUMOUR')
    normal_geno = parse_genotypes(vcf_line, method, 'NORMAL')
    info = {
        'ref': vcf_line['REF'],
        'alt': vcf_line['ALT'],
        'filter': vcf_line.get('FILTER', None),
        'methods': [method.upper()],
        't_depth': tumor_geno.depth_of_coverage,
        't_ref_count': tumor_geno.ref_depth,
        't_alt_count': tumor_geno.alt_depth,
        'n_depth': normal_geno.depth_of_coverage,
        'n_ref_count': normal_geno.ref_depth,
        'n_alt_count': normal_geno.alt_depth
    }
    return info


def transform(vcf_dir="source/gdsc/vcfs/*",
              cellline_lookup_path="source/ccle/cellline_id_lookup.tsv",
              emitter_name="json",
              emitter_prefix=None,
              emitter_directory="../pre-outputs/gdsc"):

    id_lookup = bmeg.ioutils.read_lookup(cellline_lookup_path)

    emitter = new_emitter(name=emitter_name, directory=emitter_directory, prefix=emitter_prefix)

    vcf_files = [f for f in glob(vcf_dir) if os.path.isfile(f)]
    vcf_files.sort()
    logging.debug('vcf files: {}'.format(vcf_files))
    total = len(vcf_files)
    if total == 0:
        raise ValueError("no VCF files found for: {}".format(vcf_dir))

    emitted_alleles = {}
    emitted_callsets = {}
    i = 1
    for f in vcf_files:
        logging.info('processing vcf {}/{}: {}'.format(i, total, f))
        sample_name = os.path.basename(f).split('.')[0]
        method = os.path.basename(f).split('.')[1]
        tumor_aliquot_id = id_lookup.get(sample_name, sample_name)
        vcf = read_vcf(f)
        for index, line in vcf.iterrows():
            line = line.dropna().to_dict()

            allele = make_minimal_allele(line)
            if allele.gid() not in emitted_alleles:
                emitter.emit_vertex(allele)
                emitted_alleles[allele.gid()] = None

            callset = SomaticCallset(
                id=SomaticCallset.make_gid(PROJECT_NAME, tumor_aliquot_id, None),
                tumor_aliquot_id=tumor_aliquot_id,
                normal_aliquot_id=None,
                project_id=Project.make_gid(PROJECT_NAME)
            )
            if callset.gid() not in emitted_callsets:
                emitter.emit_vertex(callset)
                emitted_callsets[callset.gid()] = None
                emitter.emit_edge(
                    SomaticCallset_Aliquots_Aliquot(
                        from_gid=callset.gid(),
                        to_gid=Aliquot.make_gid("%s:%s" % (PROJECT_NAME, callset.tumor_aliquot_id)),
                    ),
                    emit_backref=True
                )

            emitter.emit_edge(
                SomaticCallset_Alleles_Allele(
                    from_gid=callset.gid(),
                    to_gid=allele.gid(),
                    data=make_variant_call_data(line, method=method)
                ),
                emit_backref=True
            )
        i += 1

    emitter.close()


if __name__ == '__main__':  # pragma: no cover
    parser = default_argument_parser(emitter_directory_default='../pre-outputs/gdsc',
                                     emitter_prefix_default=None)
    parser.add_argument('--vcf-pattern', type=str,
                        help='vertex file pattern to glob for [e.g. source/gdsc/vcfs/*.pindel.annot.vcf.gz]',
                        required=True)
    options = parser.parse_args()
    default_logging(options.loglevel)

    if not options.emitter_prefix:
        raise ValueError("must provide value for --emitter-prefix")

    transform(vcf_dir=options.vcf_pattern,
              emitter_directory=options.emitter_directory,
              emitter_prefix=options.emitter_prefix,
              emitter_name=options.emitter)
