""" drug name, return Compound """

from bmeg.vertex import Compound


def enrich(name):
    """ retrieve payload from myvariant.info location query"""
    # TODO - this is just a stub for now
    return Compound(term_id='TODO~{}'.format(name), term=name)
