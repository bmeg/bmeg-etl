
import hashlib

class GID(str):
    def __new__(cls, content):
        return str.__new__(cls, content)

def make_allele_gid(genome, chromosome, start, end, ref_bases, alt_bases):
    # TODO
    # figure out better hashing strategy
    vid = "%s:%s:%d:%d:%s:%s" % (genome, chromosome,
                                 start, end, ref_bases,
                                 alt_bases)
    vid = vid.encode('utf-8')
    vidhash = hashlib.sha1()
    vidhash.update(vid)
    vidhash = vidhash.hexdigest()
    return "%s:%s" % ("Allele", vidhash)
