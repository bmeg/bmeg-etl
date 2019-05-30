import argparse
import logging


def default_args():
    logging.basicConfig(level=logging.INFO)
    parser = argparse.ArgumentParser()
    parser.add_argument("--project_id")
    parser.add_argument("--vertex_name")
    args = parser.parse_args()
    return args
