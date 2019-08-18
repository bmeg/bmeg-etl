from transform.gdc.gdc_transformer import GDCTransformer
from collections import defaultdict
from collections import defaultdict

EXPECTED = [('Aliquot', 182042),('Case', 33605), ('Demographic', 33340), ('Diagnosis', 33371),('Project', 47), ('Sample', 60137), ('Source', 1)]


def test_transformer():
    gdc_transformer = GDCTransformer()
    gdc = gdc_transformer.graph()
    
    lable_counts = defaultdict(int)

    for k,v in gdc.nodes.data():
        lable_counts[v['label']] += 1

    actual = [(k, lable_counts[k]) for k in sorted(lable_counts)]
    assert actual == EXPECTED, 'Should have expected label counts'

