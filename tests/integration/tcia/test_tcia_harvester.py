from transform.tcia.tcia_harvester import load_all
# from transform.tcia.tcia_harvester import harvest
# import os


def test_harvester():
    """TODO harvest gets all the data, should mock it for test"""
    # harvest(os.environ['TCIA_APIKEY'])
    (all_patients, all_patient_study, all_series) = load_all()
    for a in [all_patients, all_patient_study, all_series]:
        assert len(a) > 0, 'Should be > 0'
