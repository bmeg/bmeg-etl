
""" utility, minimal logging   """
import logging
import json


def default_logging(loglevel):
    """ bare bones logging"""
    fmt = '%(asctime)s - %(levelname)s - %(threadName)s - %(message)s'
    logging.basicConfig(level=loglevel, format=fmt)


def log_missing_vertex(msg):
    """ standard way to log missing data"""
    msg['DATA_EXCEPTION'] = 'MISSING_VERTEX'
    logging.warning(json.dumps(msg, separators=(',', ':')))
