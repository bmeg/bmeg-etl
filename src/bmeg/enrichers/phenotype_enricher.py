""" given disease name, return Phenotype """

from bmeg.vertex import Phenotype


def phenotype_factory(name):
    """ create a stub compound for downstream normalization """
    return Phenotype(term_id='TODO:{}'.format(name), term='TODO', name=name)
