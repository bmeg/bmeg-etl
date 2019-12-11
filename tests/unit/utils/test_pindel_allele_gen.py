import os
import pandas
import pytest

from bmeg.vcf import make_minimal_pindel_allele
from bmeg.maf import make_minimal_allele


@pytest.fixture
def min_maf_path(request):
    """unnormalized minimal maf"""
    return os.path.join(request.fspath.dirname, 'min.maf')


@pytest.fixture
def anno_maf_path(request):
    """normalized maf"""
    return os.path.join(request.fspath.dirname, 'anno.maf')


def test_make_minimal_pindel_allele(min_maf_path, anno_maf_path):
    min_maf = pandas.read_csv(min_maf_path, sep='\t', comment='#', dtype=str)
    anno_maf = pandas.read_csv(anno_maf_path, sep='\t', comment='#', dtype=str)
    for index, line in min_maf.iterrows():
        min_line = line.dropna().to_dict()
        min_vcf = {'#CHROM': min_line["Chromosome"],
                   "POS": min_line["Start_Position"],
                   "REF": min_line["Reference_Allele"],
                   "ALT": min_line["Tumor_Seq_Allele2"]}
        anno_line = anno_maf.iloc[index, :].dropna().to_dict()
        anno_allele = make_minimal_allele(anno_line)
        min_allele = make_minimal_pindel_allele(min_vcf)

        if min_allele.gid() != anno_allele.gid():
            print("original:", min_allele.gid(), min_allele.props())
            print("normalized:", anno_allele.gid(), anno_allele.props())
            raise AssertionError("gids do not match")
