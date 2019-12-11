import pandas

from bmeg import Allele, Project
from bmeg.ioutils import reader
chromosome = '#CHROM'
start = 'POS'
reference_bases = 'REF'
alternate_bases = 'ALT'


def read_vcf(filename):
    comments = 0
    with reader(filename) as fh:
        for line in fh:
            if line.startswith('#CHROM'):
                break
            else:
                comments += 1
    return pandas.read_table(filename, skiprows=comments, header=0, dtype=str)


def make_minimal_allele(vcf_line, genome='GRCh37'):
    assert isinstance(vcf_line, dict)
    if reference_bases not in vcf_line and alternate_bases not in vcf_line:
        raise ValueError('ref and alt bases are missing for row: %s', vcf_line)
    if reference_bases not in vcf_line:
        vcf_line[reference_bases] = '.'
    if alternate_bases not in vcf_line:
        vcf_line[alternate_bases] = '.'
    return Allele(
        id=Allele.make_gid(
            genome, vcf_line.get(chromosome), vcf_line.get(start),
            vcf_line.get(reference_bases), vcf_line.get(alternate_bases),
        ),
        chromosome=vcf_line.get(chromosome),
        start=int(vcf_line.get(start)),
        reference_bases=vcf_line.get(reference_bases),
        alternate_bases=vcf_line.get(alternate_bases),
        genome=genome,
        project_id=Project.make_gid('Reference')
    )


def make_minimal_pindel_allele(vcf_line, genome='GRCh37'):
    allele = make_minimal_allele(vcf_line, genome)

    pos = allele.start
    ref = allele.reference_bases
    alt = allele.alternate_bases

    if len(ref) > len(alt):
        ref = ref[1:]
        if len(alt) == 1:
            alt = "-"
        else:
            alt = alt[1:]
        pos += 1
    elif len(ref) < len(alt):
        alt = alt[1:]
        if len(ref) == 1:
            ref = "-"
        else:
            ref = ref[1:]

    allele.id = Allele.make_gid(allele.genome, allele.chromosome, pos, ref, alt)

    return allele
