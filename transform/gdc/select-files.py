import os
import sys
import json
import argparse

from sets import Set
from shutil import copyfile

def find_files(names, path):
    name_set = Set(names)
    found = []

    for root, dirs, files in os.walk(path):
        file_set = Set(files)
        intersect = name_set & file_set
        absolutes = [os.path.join(root, f) for f in intersect]
        found = found + absolutes

    return found

def select_files(options):
    with open(options.files) as j:
        files = json.loads(j.read())

    found = find_files(files, options.path)
    for file in found:
        copyfile(file, os.path.join(options.target, os.path.basename(file)))

def parse_args(args):
    args = args[1:]
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('--files', type=str, help='a json file describing the files to select')
    parser.add_argument('--path', type=str, help='path to a directory tree containing the files')
    parser.add_argument('--target', type=str, help='path to the directory to store the selected files')

    return parser.parse_args(args)

if __name__ == '__main__':
    options = parse_args(sys.argv)
    select_files(options)