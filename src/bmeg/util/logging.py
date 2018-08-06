#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" utility, minimal logging   """
import logging


def default_logging(loglevel):  # pragma: no cover
    """ bear bones logging"""
    fmt = '%(asctime)s - %(levelname)s - %(threadName)s - %(message)s'
    logging.basicConfig(level=loglevel, format=fmt)
