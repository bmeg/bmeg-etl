from __future__ import print_function

import argparse
import os
import ujson
from types import SimpleNamespace as SN

import bmeg.ioutils


def summarize_output(path, emit_ids):
    """summarize missing vertices by label"""
    if not os.path.isfile(path):
        raise ValueError("path {} is not a file".format(path))

    labels = {}
    input_stream = bmeg.ioutils.reader(path)
    for line in input_stream:
        err = SN(**ujson.loads(line))
        elabel = err.Label
        if err.msg == "From does not exist":
            label = err.From.split(":")[0]
            gid = err.From
        elif err.msg == "To does not exist":
            label = err.To.split(":")[0]
            gid = err.To
        else:
            raise ValueError("unknown error message: {}".format(err.msg))

        if label.startswith("ENSP"):
            label = "Protein"
        elif label.startswith("ENST"):
            label = "Transcript"
        elif label.startswith("ENSG"):
            label = "Gene"

        if label not in labels:
            labels[label] = {}
            labels[label][elabel] = [gid]
        else:
            if elabel not in labels[label]:
                labels[label][elabel] = [gid]
            else:
                if gid not in labels[label][elabel]:
                    labels[label][elabel].append(gid)

    if emit_ids:
        print(ujson.dumps(labels))
        return

    labelCounts = {}
    for k, v in labels.items():
        labelCounts[k] = {}
        for ek, ev in v.items():
            labelCounts[k][ek] = len(ev)
    print(ujson.dumps(labelCounts))
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("path",
                        help="path to check graph output")
    parser.add_argument("--emit-gids",
                        default=False,
                        action="store_true",
                        help="emit missing gids rather than label counts")
    args = parser.parse_args()
    summarize_output(args.path, args.emit_gids)
