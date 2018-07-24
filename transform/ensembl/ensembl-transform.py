import argparse
import csv
import gzip
import io
import os
import re

from urllib.parse import unquote

from bmeg.models.vertex_models import Gene, Transcript, Exon
from bmeg.models.edge_models import TranscriptFor, ExonFor
from bmeg.models.emitter import Emitter


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
    emitter = Emitter(args.output_prefix)
    emit = emitter.emit

    if args.gz:
        inhandle = io.TextIOWrapper(gzip.GzipFile(args.input))
    else:
        inhandle = open(args.input, "r")

    reader = csv.DictReader(
        filter(lambda row: row[0] != '#', inhandle),
        delimiter="\t",
        fieldnames=["seqId", "source", "type", "start", "end", "score",
                    "strand", "phase", "attributes"]
    )

    for line in reader:
        attrs = parse_attributes(line["attributes"])

        if line["type"] == "gene":
            g = Gene(ensembl_id=attrs["gene_id"],
                     symbol=attrs["Name"],
                     description=attrs.get("description", ""),
                     chromosome=line["seqId"],
                     start=int(line["start"]),
                     end=int(line["end"]),
                     strand=line["strand"],
                     mygeneinfo=None)
            emit(g)

        elif line["type"] == "transcript" or line["type"] == "mRNA":
            t = Transcript(ensembl_id=attrs["transcript_id"],
                           gene_id=get_parent_gene(attrs["Parent"]),
                           chromosome=line["seqId"],
                           start=int(line["start"]),
                           end=int(line["end"]),
                           strand=line["strand"])
            emit(t)

        elif line["type"] == "exon":
            e = Exon(ensembl_id=attrs["exon_id"],
                     transcript_id=get_parent_transcript(attrs["Parent"]),
                     chromosome=line["seqId"],
                     start=int(line["start"]),
                     end=int(line["end"]),
                     strand=line["strand"])
            emit(e)

    emitter.close()
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', "-i", type=str, required=True,
                        help='Path to the Ensembl GFF3 file')
    parser.add_argument('--gz', action="store_true", default=False,
                        help='Input file is gzipped')
    parser.add_argument('--output-prefix', "-o", type=str, required=True,
                        help='Output file prefix')
    args = parser.parse_args()
    args.output_prefix = os.path.abspath(args.output_prefix)
    transform(args)
