from array import array
import csv
import os
import struct


def default_gene_func(name, description):
    return name


def split_ensembl_id(name, description):
    # Assumes the "name" is an ensembl gene ID, such as "ENSG00000240361.1".
    return name.split(".")[0]


def parse_gct(path, output_dir, gene_func=default_gene_func):
    """
    parse_gct parses a GCT file (expression data) and yields a dictionary
    of gene expression values for each sample.

    gene_func can be used to transform the gene ID.

    https://software.broadinstitute.org/cancer/software/genepattern/file-formats-guide#GCT
    """

    # Expression values are 32-bit floats.
    cell_type = "f"
    # Expression values take 4 bytes.
    cell_bytesize = 4

    reader = csv.reader(open(path), delimiter="\t")

    # Skip the version and shape headers.
    next(reader)
    next(reader)

    # Get the samples from the column names.
    # The first two columns are "Name" and "Description", skip those.
    samples = next(reader)[2:]

    genes = []

    # GCT files can be massive, so they take up a lot of memory (many gigabytes)
    # and they're slow to access for BMEG. To avoid memory usage, the code below
    # writes a temporary, binary-formatted file that fits BMEG's access pattern.

    # write values out to temp matrix
    t = open(os.path.join(output_dir, "tmp_matrix.data"), "wb")
    for row in reader:
        gene = row[0]
        desc = row[1]
        genes.append(gene_func(gene, desc))

        values = row[2:]
        assert(len(samples) == len(values))

        v = array(cell_type, (float(a) for a in values))
        v.tofile(t)
    t.close()

    # read back, steping through the rows to do a transpose
    t = open(os.path.join(output_dir, "tmp_matrix.data"), "rb")
    rowSize = cell_bytesize * len(samples)
    for i, sample in enumerate(samples):
        values = {}
        for j, gene in enumerate(genes):
            t.seek(i * cell_bytesize + j * rowSize)
            o = struct.unpack(cell_type, t.read(cell_bytesize))
            values[gene] = o[0]

        yield sample, values

    t.close()
    os.unlink(os.path.join(output_dir, "tmp_matrix.data"))
