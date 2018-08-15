
""" transform a maf file into vertexs[variant, allele]   """
from bmeg.vertex import Biosample, Callset, Gene
from bmeg.edge import AlleleCall

from bmeg.maf.maf_transform import main, get_value,  MAFTransformer
from bmeg.maf.maf_transform import transform as parent_transform

TUMOR_SAMPLE_BARCODE = "Tumor_Sample_Barcode"  # 15
NORMAL_SAMPLE_BARCODE = "Matched_Norm_Sample_Barcode"  # 16


class MC3_MAFTransformer(MAFTransformer):

    def barcode_to_sampleid(self, barcode):
        """ create tcga sample barcode """
        return "-".join(barcode.split("-")[0:4])

    def create_gene_gid(self, line):  # pragma nocover
        """ override, create gene_gid from line """
        ensembl_id = line.get('Gene', None)
        return Gene.make_gid(gene_id=ensembl_id)

    def allele_call_maker(self, allele, line=None):
        """ create call from line """
        keys = ['t_depth', 't_ref_count', 't_alt_count', 'n_depth',
                'n_ref_count', 'n_alt_count', 'FILTER',
                'Match_Norm_Seq_Allele1',	'Match_Norm_Seq_Allele2',
                'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2',
                ]
        info = {}
        for k in keys:
            info[k] = get_value(line, k, None)
        return AlleleCall(info)

    def callset_maker(self, allele, source, centerCol, method, line):
        """ create callset from line """
        barcode = line[TUMOR_SAMPLE_BARCODE]
        sample = barcode
        sample = self.barcode_to_sampleid(barcode)
        sample = Biosample.make_gid(sample)
        sample_callsets = []
        sample_calls = []
        if centerCol in line:
            for c in line[centerCol].split("|"):
                center = c.replace("*", "")
                # callset_id = "%s:%s" % (sample, center)
                callset = Callset(sample,
                                  line[NORMAL_SAMPLE_BARCODE],
                                  center, source)
                sample_callsets.append(callset)
                sample_calls.append((self.allele_call_maker(allele, line),
                                     callset.gid()))
        else:
            callset = Callset(sample,
                              line[NORMAL_SAMPLE_BARCODE],
                              method, source)
            sample_callsets.append(callset)
            sample_calls.append((self.allele_call_maker(allele, line),
                                 callset.gid()))
        return sample_calls, sample_callsets


def transform(mafpath, prefix, workers=5, skip=0, harvest=True, filter=[]):
    return parent_transform(mafpath, prefix, workers, skip, harvest, filter,
                            transformer=MC3_MAFTransformer())


if __name__ == '__main__':  # pragma: no cover
    main()
