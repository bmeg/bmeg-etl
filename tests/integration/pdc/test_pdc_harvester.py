from transform.pdc.pdc_harvester import load_all


def test_harvester():
    """TODO harvest gets all the data, should mock it for test"""
    # harvest()
    (all_cases, all_programs, all_studies) = load_all()
    for a in [all_cases, all_programs, all_studies]:
        assert len(a) > 0, 'Should be > 0'
