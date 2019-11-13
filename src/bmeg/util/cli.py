""" utility, sets common arguments   """
import argparse
import logging


def default_argument_parser(emitter_directory_default=None, emitter_prefix_default=None):  # pragma: no cover
    # We don't need the first argument, which is the program name
    # Construct the parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        '--emitter-prefix', type=str,
        default=emitter_prefix_default,
        help='filename prefix for emitter output files'
    )
    parser.add_argument(
        '--emitter-directory', type=str,
        default=emitter_directory_default,
        help='emitter output directory'
    )
    parser.add_argument(
        '-d', '--debug',
        help="Print lots of debugging statements",
        action="store_const", dest="loglevel", const=logging.DEBUG,
        default=logging.WARNING,
    )
    parser.add_argument(
        '-v', '--verbose',
        help="Be verbose",
        action="store_const", dest="loglevel", const=logging.INFO,
    )
    parser.add_argument('--emitter', type=str,
                        default='json',
                        choices=["json", "debug"],
                        help='emitter type [json, debug]')
    return parser
