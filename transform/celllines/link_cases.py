import json

import bmeg.ioutils
from bmeg.emitter import JSONEmitter
from bmeg import (Case,
                  Case_SameAs_Case)


def transform(case_files=["outputs/ccle/ccle.Case.Vertex.json.gz",
                          "outputs/ctrp/ctrp.Case.Vertex.json.gz",
                          "outputs/gdsc/gdsc.Case.Vertex.json.gz"],
              emitter_prefix=None,
              emitter_directory="celllines"):

    emitter = JSONEmitter(directory=emitter_directory, prefix=emitter_prefix)
    emitted_edges = {}
    cases = {}
    for path in case_files:
        with bmeg.ioutils.reader(path) as f:
            for line in f:
                case = json.loads(line)
                case_id = case["data"]["case_id"]
                if case_id not in cases:
                    cases[case_id] = [case["gid"]]
                else:
                    cases[case_id] = cases[case_id] + [case["gid"]]

    for k, v in cases.items():
        i = 0
        while i < len(v) - 1:
            j = i + 1
            f = v[i]
            t = v[j]
            if (f, t) not in emitted_edges:
                emitter.emit_edge(
                    Case_SameAs_Case(
                        from_gid=Case._gid_cls(f),
                        to_gid=Case._gid_cls(t)
                    ),
                )
                emitted_edges[(f, t)] = True
            if (t, f) not in emitted_edges:
                emitter.emit_edge(
                    Case_SameAs_Case(
                        from_gid=Case._gid_cls(t),
                        to_gid=Case._gid_cls(f)
                    ),
                )
                emitted_edges[(t, f)] = True
            i += 1
    emitter.close()
    return


if __name__ == '__main__':  # pragma: no cover
    transform()
