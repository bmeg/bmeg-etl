#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" utility, minimal logging   """
from logging import *


def default_logging(loglevel):  # pragma: no cover
    """ bear bones logging"""
    fmt = '%(asctime)s - %(levelname)s - %(threadName)s - %(message)s'
    basicConfig(level=loglevel, format=fmt)


def simple_logging():
    fmt = '%(levelname)s %(message)s'
    basicConfig(level=INFO, format=fmt)
