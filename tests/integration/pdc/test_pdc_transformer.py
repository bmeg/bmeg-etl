from transform.pdc.pdc_transformer import PDCTransformer
from collections import defaultdict


def test_transformer():
    pdc_transformer = PDCTransformer()

    for a in [pdc_transformer.all_cases, pdc_transformer.all_programs, pdc_transformer.all_studies]:
        assert len(a) > 0, 'Should be > 0'

    assert len(pdc_transformer.sample_gdc_projects.keys()) == 4, 'Should have 4 projects'

    G = pdc_transformer.graph()
    assert len(G.nodes) == 42472, 'Should have lots of nodes'

    node_counts = defaultdict(int)
    for k, v in list(G.nodes.data()):
        node_counts[v['label']] += 1
    actual = sorted(node_counts.items())
    print(actual)
    expected = [('Aliquot', 1725), ('Case', 977), ('Demographic', 977), ('Diagnosis', 977), ('File', 36461), ('Project',4), ('Sample', 1330), ('Source', 1), ('Study', 20)]
    assert actual == expected, 'Graph should have expected number of nodes'
