def set_gid(obj, gid):
    object.__setattr__(obj, "gid", gid)


def get_tcga_individual_barcode(id):
    parts = id.split("-")
    return "-".join(parts[0:3])


def get_tcga_sample_barcode(id):
    parts = id.split("-")
    return "-".join(parts[0:4])


def get_tcga_portion_barcode(id):
    parts = id.split("-")
    parts[5] = parts[5][:-1]
    return "-".join(parts[0:5])


def get_tcga_analyte_barcode(id):
    parts = id.split("-")
    return "-".join(parts[0:5])


def get_tcga_aliquot_barcode(id):
    parts = id.split("-")
    return "-".join(parts[0:7])


def tcga_barcode_is_tumor(id):
    parts = id.split("-")
    sample_number = parts[4][:-1]
    return sample_number < 10


def tcga_barcode_is_normal(id):
    parts = id.split("-")
    sample_number = parts[4][:-1]
    return sample_number >= 10
