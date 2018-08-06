#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" given a genomic location, harvest allele and dependencies """
import requests
import urllib
import logging
import json
from bmeg.vertex import Allele


TIMEOUT = 30


def _myvariantinfo(genome, chromosome, start, end, reference_bases,
                   alternate_bases, annotations=[]):
    """ retrieve payload from myvariant.info location query"""
    if (reference_bases == '' or alternate_bases == '' or
       reference_bases is None or alternate_bases is None):
            raise ValueError('reference_bases & alternate_bases must be set')
    url = hit = None
    try:
        dbSNP = None
        for name in annotations:
            # "dbSNP_RS:novel"
            if 'dbSNP_RS' not in name:
                continue
            if name == 'dbSNP_RS:novel':
                continue
            if name == 'dbSNP_RS:.':
                continue
            dbSNP = name.split(':')[1]

        vcf_alternate = alternate_bases
        if alternate_bases == '-':
            vcf_alternate = reference_bases

        q = 'chrom:{} AND hg19.start:{} AND hg19.end:{} AND vcf.alt:{}' \
            .format(chromosome, start, end, vcf_alternate)
        # logging.debug(genome, chromosome, start, end, alternate_bases)
        logging.debug(q)
        url = 'http://myvariant.info/v1/query?q={}' \
            .format(urllib.parse.quote_plus(q))
        response = requests.get(url, timeout=TIMEOUT)
        logging.debug(response)
        _status_code = response.status_code
        if _status_code == 200:
            r = response.json()
            logging.debug('hit length = {}'.format(len(r['hits'])))
            # return first response
            for hit in r['hits']:
                return hit

        # no hit?, lookup by snp if we have it
        # ["00058da0fb86f362e1504b6ec02d5e4446d9db44","GRCh37","16","397035","397035","C","T",["Variant_Type:SNP","Feature_type:Transcript","Feature:ENST00000262320","dbSNP_RS:rs375097914"]]
        if dbSNP:
            url = 'http://myvariant.info/v1/variant/{}'.format(dbSNP)
            response = requests.get(url, timeout=TIMEOUT)
            _status_code = response.status_code
            if _status_code == 200:
                hits = response.json()
                if not type(hits) is list:
                    hits = [hits]
                for hit in hits:
                    # ensure we  have an exact match
                    if (hit['vcf']['alt'] == alternate_bases and
                        (hit['vcf']['ref'] == '-' or
                            hit['vcf']['ref'] == reference_bases)):
                        # TODO - can this happen? do we have a test case?
                        return hit
                    else:
                        _log_json({
                                        'stage': 'dbSNP_mismatch',
                                        'genome': genome,
                                        'chromosome': chromosome,
                                        'start': start,
                                        'end': end,
                                        'reference_bases': reference_bases,
                                        'alternate_bases': alternate_bases,
                                        'annotations': annotations,
                                        'url': url,
                                        'hit': hit
                                    })
    except Exception as e:  # pragma: no cover
        logging.exception(json.dumps({
                        'stage': 'myvariantinfo_exception',
                        'exception': str(e),
                        'genome': genome,
                        'chromosome': chromosome,
                        'start': start,
                        'end': end,
                        'reference_bases': reference_bases,
                        'alternate_bases': alternate_bases,
                        'annotations': annotations,
                        'url': url,
                        'hit': hit
                    }))
        return None
    # default
    _log_json({
                    'stage': 'myvariantinfo_nofind',
                    'genome': genome,
                    'chromosome': chromosome,
                    'start': start,
                    'end': end,
                    'reference_bases': reference_bases,
                    'alternate_bases': alternate_bases,
                    'annotations': annotations,
                })


def harvest(genome, chromosome, start, end, reference_bases, alternate_bases,
            annotations=[], harvest=True, filter={}):
    """ creates an Allele, including external datasources (myvariant.info)"""
    allele = None
    allele = Allele(genome, chromosome, start, end, reference_bases,
                    alternate_bases, annotations)
    if allele.gid() in filter:
        _log_json({
                        'stage': 'filtered',
                        'allele_gid':  allele.gid(),
                        'genome': genome,
                        'chromosome': chromosome,
                        'start': start,
                        'end': end,
                        'reference_bases': reference_bases,
                        'alternate_bases': alternate_bases,
                        'annotations': annotations,
                    })
        return None

    if harvest:
        myvariantinfo_dict = _myvariantinfo(genome, chromosome, start, end,
                                            reference_bases, alternate_bases,
                                            annotations)
        if myvariantinfo_dict:
            allele = Allele(genome, chromosome, start, end, reference_bases,
                            alternate_bases, annotations, myvariantinfo_dict)
    return allele


def _log_json(obj):
    ''' log object as json '''
    logging.info(json.dumps(obj))


# def create(genome, chromosome, start, end, reference_bases, alternate_bases,
#            annotations):
#     """ creates an Allele, no external datasources"""
#     return Allele(genome, chromosome, start, end, reference_bases,
#                   alternate_bases, annotations)
