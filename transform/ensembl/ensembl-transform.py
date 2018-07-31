import argparse
import csv
import os
import re

from urllib.parse import unquote

import bmeg.ioutils
from bmeg.vertex import Gene, Transcript, Exon
from bmeg.edge import TranscriptFor, ExonFor
from bmeg.emitter import JSONEmitter


def parse_attributes(attrs):
    parsed = {}
    attrs = attrs.split(";")
    for a in attrs:
        a = a.split("=")
        key = a[0]
        val = a[1].split(",")
        if len(val) > 1:
            parsed[key] = [unquote(x) for x in val]
        else:
            parsed[key] = unquote(val[0])
    return parsed


def get_parent_gene(parent):
    if isinstance(parent, str):
        p = re.sub("^gene:", "", parent)
        if not p.startswith("ENSG"):
            raise ValueError("%s does not look like a gene_id" % parent)
        return p
    else:
        raise TypeError("expected 'Parent' to be a str")


def get_parent_transcript(parent):
    if isinstance(parent, str):
        p = re.sub("^transcript:", "", parent)
        if not p.startswith("ENST"):
            raise ValueError("%s does not look like a transcript_id" % parent)
        return p
    else:
        raise TypeError("%s is not a str" % parent)


def transform(args):
    """
    Transform the file downloaded from:
        ftp://ftp.ensembl.org/pub/grch37/update/gff3/homo_sapiens/Homo_sapiens.GRCh37.87.gff3.gz
    """
    emitter = JSONEmitter(args.output_prefix)
    emit = emitter.emit

    inhandle = bmeg.ioutils.reader(args.input)
    reader = csv.DictReader(
        filter(lambda row: row[0] != '#', inhandle),
        delimiter="\t",
        fieldnames=["seqId", "source", "type", "start", "end", "score",
                    "strand", "phase", "attributes"]
    )

    for line in reader:
        attrs = parse_attributes(line["attributes"])

        if line["type"] == "gene":
            g = Gene(gene_id=attrs["gene_id"],
                     symbol=attrs["Name"],
                     description=attrs.get("description", ""),
                     chromosome=line["seqId"],
                     start=int(line["start"]),
                     end=int(line["end"]),
                     strand=line["strand"],
                     mygeneinfo=None)
            emit(g)

        elif line["type"] == "transcript" or line["type"] == "mRNA":
            gene_id = get_parent_gene(attrs["Parent"])
            t = Transcript(transcript_id=attrs["transcript_id"],
                           gene_id=gene_id,
                           chromosome=line["seqId"],
                           start=int(line["start"]),
                           end=int(line["end"]),
                           strand=line["strand"])
            emit(t)

            tf = TranscriptFor(from_gid=t.gid,
                               to_gid=Gene.make_gid(gene_id))
            emit(tf)

        elif line["type"] == "exon":
            transcript_id = get_parent_transcript(attrs["Parent"])
            e = Exon(exon_id=attrs["exon_id"],
                     transcript_id=transcript_id,
                     chromosome=line["seqId"],
                     start=int(line["start"]),
                     end=int(line["end"]),
                     strand=line["strand"])
            emit(e)

            ef = ExonFor(from_gid=e.gid,
                         to_gid=Transcript.make_gid(transcript_id))
            emit(ef)

    emitter.close()
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', "-i", type=str, required=True,
                        help='Path to the Ensembl GFF3 file')
    parser.add_argument('--output-prefix', "-o", type=str, required=True,
                        help='Output file prefix')
    args = parser.parse_args()
    args.output_prefix = os.path.abspath(args.output_prefix)
    transform(args)
