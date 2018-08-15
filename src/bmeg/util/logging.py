#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" utility, minimal logging   """
import logging


# shortcuts for code importing this module
error = logging.error
warn = logging.warn
info = logging.info
debug = logging.debug


def default_logging(loglevel):  # pragma: no cover
    """ bare bones logging"""
    fmt = '%(asctime)s - %(levelname)s - %(threadName)s - %(message)s'
    logging.basicConfig(level=loglevel, format=fmt)


def simple_logging():
    fmt = '%(levelname)s %(message)s'
    logging.basicConfig(level=logging.INFO, format=fmt)
