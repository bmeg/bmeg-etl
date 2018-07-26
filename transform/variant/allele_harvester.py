#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" given a genomic location, harvest allele and dependencies """
import requests
import urllib
import logging

from bmeg.models.vertex_models import Allele


def _myvariantinfo(genome, chromosome, start, end, reference_bases,
                   alternate_bases, annotations):
    """ retrieve payload from myvariant.info location query"""
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
        url = 'http://myvariant.info/v1/query?q={}' \
            .format(urllib.parse.quote_plus(q))
        response = requests.get(url)
        _status_code = response.status_code
        if _status_code == 200:
            r = response.json()
            # return first response
            for hit in r['hits']:
                return hit

        # no hit?, lookup by snp if we have it
        # ["00058da0fb86f362e1504b6ec02d5e4446d9db44","GRCh37","16","397035","397035","C","T",["Variant_Type:SNP","Feature_type:Transcript","Feature:ENST00000262320","dbSNP_RS:rs375097914"]]
        if dbSNP:
            url = 'http://myvariant.info/v1/variant/{}'.format(dbSNP)
            response = requests.get(url)
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
                        return hit
                    else:
                        logging.info({
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
    except Exception as e:
        logging.exception({
                        'stage': 'myvariantinfo_exception',
                        'message': e.message,
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
        return None
    # default
    logging.debug({
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
            annotations):
    """ creates an Allele, including external datasources (myvariant.info)"""
    myvariantinfo_dict = _myvariantinfo(genome, chromosome, start, end,
                                        reference_bases, alternate_bases,
                                        annotations)
    if myvariantinfo_dict:
        return Allele(genome, chromosome, start, end, reference_bases,
                      alternate_bases, annotations, myvariantinfo_dict)
    return Allele(genome, chromosome, start, end, reference_bases,
                  alternate_bases, annotations)
