#!/usr/bin/env python

import gzip
import json
import sys
import torch
from scipy.spatial.distance import cdist

useTorch = True

def load_file(path):
    data = {}
    with gzip.open(path) as handle:
        for line in handle:
            m = json.loads(line)
            gid = m["gid"]
            if "morgan_fingerprint_2" in m["data"]:
                n = m["data"]["morgan_fingerprint_2"]
                data[gid] = n
    return data


if __name__ == "__main__":
    data = load_file(sys.argv[1])
    if useTorch:
        X = torch.tensor(list(data.values()), dtype=torch.float32)
        y = torch.cdist(X, X)
    else:
        X = list(data.values())
        y = cdist(X,X)
    print(y)
