import csv
import re

from urllib.parse import unquote

from bmeg import Exon, Gene, Transcript

import bmeg.ioutils
from bmeg.emitter import JSONEmitter


GENOME_BUILD = "GRCh37"
DEFAULT_DIRECTORY = "ensembl"


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


def transform(
    emitter_directory=DEFAULT_DIRECTORY,
    gff3_path="source/ensembl/Homo_sapiens.GRCh37.87.chr_patch_hapl_scaff.gff3.gz"
):
    """
    Transform the file downloaded from:
        ftp://ftp.ensembl.org/pub/grch37/release-94/gff3/homo_sapiens/Homo_sapiens.GRCh37.87.chr_patch_hapl_scaff.gff3.gz
    """

    emitter = JSONEmitter(directory=emitter_directory, prefix=None)

    inhandle = bmeg.ioutils.reader(gff3_path)
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

        if attrs["feature_id"] is not None:
            if attrs["feature_id"] not in features:
                features[attrs["feature_id"]] = [attrs]
            else:
                features[attrs["feature_id"]].append(attrs)

    emitted_genes = {}
    emitted_transcripts = {}
    for key, aset in features.items():
        if aset[0]["type"] == "exon":
            transcripts = []
            for attrs in aset:
                transcript_id = get_parent_transcript(attrs["Parent"])
                transcripts.append(transcript_id)
            attrs = aset[0]
            e = Exon(exon_id=attrs["exon_id"],
                     transcript_id=transcripts,
                     chromosome=attrs["seqId"],
                     start=attrs["start"],
                     end=attrs["end"],
                     strand=attrs["strand"],
                     genome=GENOME_BUILD)
            emitter.emit(e)

            for transcript_id in transcripts:
                if transcript_id not in emitted_transcripts:
                    for attrs in features[transcript_id]:
                        gene_id = get_parent_gene(attrs["Parent"])
                        t = Transcript(transcript_id=attrs["transcript_id"],
                                       gene_id=gene_id,
                                       chromosome=attrs["seqId"],
                                       start=int(attrs["start"]),
                                       end=int(attrs["end"]),
                                       strand=attrs["strand"],
                                       biotype=attrs["type"],
                                       genome=GENOME_BUILD)
                        emitter.emit(t)
                        emitted_transcripts[transcript_id] = True

                        if gene_id not in emitted_genes:
                            for attrs in features[gene_id]:
                                g = Gene(gene_id=attrs["gene_id"],
                                         symbol=attrs["Name"],
                                         description=attrs.get("description", ""),
                                         chromosome=attrs["seqId"],
                                         start=int(attrs["start"]),
                                         end=int(attrs["end"]),
                                         strand=attrs["strand"],
                                         genome=GENOME_BUILD)
                                emitter.emit(g)
                                emitted_genes[gene_id] = True
    emitter.close()


if __name__ == '__main__':  # pragma: no cover
    transform()
