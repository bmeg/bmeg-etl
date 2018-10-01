import csv
import re

from urllib.parse import unquote

import bmeg.ioutils
from bmeg.vertex import Gene, Transcript, Exon
from bmeg.edge import TranscriptFor, ExonFor
from bmeg.emitter import JSONEmitter


GENOME_BUILD = "GRCh37"


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


"""
Transform the file downloaded from:
    ftp://ftp.ensembl.org/pub/grch37/update/gff3/homo_sapiens/Homo_sapiens.GRCh37.87.gff3.gz
"""
emitter = JSONEmitter("ensembl")

inhandle = bmeg.ioutils.reader("source/ensembl/Homo_sapiens.GRCh37.87.gff3.gz")
reader = csv.DictReader(
    filter(lambda row: row[0] != '#', inhandle),
    delimiter="\t",
    fieldnames=["seqId", "source", "type", "start", "end", "score",
                "strand", "phase", "attributes"]
)

features = {}
for line in reader:
    attrs = parse_attributes(line["attributes"])
    attrs["type"] = line["type"]
    attrs["feature_id"] = attrs.get("gene_id", attrs.get("transcript_id", attrs.get("exon_id", None)))
    attrs["seqId"] = line["seqId"]
    attrs["start"] = int(line["start"])
    attrs["end"] = int(line["end"])
    attrs["strand"] = line["strand"]

    if attrs["feature_id"] is not None and attrs["feature_id"] not in features:
        features[attrs["feature_id"]] = attrs

emitted_genes = {}
emitted_transcripts = {}
for key, attrs in features.items():
    if attrs["type"] == "exon":
        transcript_id = get_parent_transcript(attrs["Parent"])
        e = Exon(exon_id=attrs["exon_id"],
                 transcript_id=transcript_id,
                 chromosome=attrs["seqId"],
                 start=attrs["start"],
                 end=attrs["end"],
                 strand=attrs["strand"],
                 genome=GENOME_BUILD)
        emitter.emit_vertex(e)
        emitter.emit_edge(ExonFor(),
                          from_gid=e.gid(),
                          to_gid=Transcript.make_gid(transcript_id))

        if transcript_id not in emitted_transcripts:
            attrs = features[transcript_id]
            gene_id = get_parent_gene(attrs["Parent"])
            t = Transcript(transcript_id=attrs["transcript_id"],
                           gene_id=gene_id,
                           chromosome=attrs["seqId"],
                           start=int(attrs["start"]),
                           end=int(attrs["end"]),
                           strand=attrs["strand"],
                           genome=GENOME_BUILD)
            emitter.emit_vertex(t)
            emitter.emit_edge(TranscriptFor(),
                              from_gid=t.gid(),
                              to_gid=Gene.make_gid(gene_id))
            emitted_transcripts[transcript_id] = True

            if gene_id not in emitted_genes:
                attrs = features[gene_id]
                g = Gene(gene_id=attrs["gene_id"],
                         symbol=attrs["Name"],
                         description=attrs.get("description", ""),
                         chromosome=attrs["seqId"],
                         start=int(attrs["start"]),
                         end=int(attrs["end"]),
                         strand=attrs["strand"],
                         genome=GENOME_BUILD)
                emitter.emit_vertex(g)
                emitted_genes[gene_id] = True

emitter.close()
