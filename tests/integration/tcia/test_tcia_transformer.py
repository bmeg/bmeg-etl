from transform.tcia.tcia_transformer import TCIATransformer
from collections import defaultdict
from collections import defaultdict

EXPECTED = [('Demographic', 18201), ('Patient', 18201), ('Project', 88), ('Series', 121158), ('Source', 1), ('Study', 31018)]


def test_transformer():
    tcia_transformer = TCIATransformer()
    tcia = tcia_transformer.graph()
    
    lable_counts = defaultdict(int)

    for k,v in tcia.nodes.data():
        lable_counts[v['label']] += 1

    actual = [(k, lable_counts[k]) for k in sorted(lable_counts)]
    print(actual)
    assert actual == EXPECTED, 'Should have expected label counts'

