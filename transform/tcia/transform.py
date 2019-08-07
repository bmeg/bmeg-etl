import json
import click
import logging
from pathlib import Path
from dotenv import find_dotenv, load_dotenv
import os
import pathlib

from tciaclient import TCIAClient


@click.command()
@click.option('-a', '--api_key', envvar='TCIA_APIKEY', default=None)
def transform(api_key):
    """ Runs data processing scripts to turn raw data from (../raw) into
        cleaned data ready to be analyzed (saved in ../processed).
    """
    assert api_key, 'Please set --api_key or export TCIA_APIKEY=...'
    logger = logging.getLogger(__name__)
    tcia_client = TCIAClient(apiKey=api_key, baseUrl="https://services.cancerimagingarchive.net/services/v3",resource = "TCIA")

    def get_patient():
        with tcia_client.get_patient() as response:
            if response.status == 200:
                data = response.read().decode('utf-8')
                patients = json.loads(data)
                c = 0
                with open('source/tcia/patients.json', 'w') as outs:
                    for p in patients:
                        json.dump(p, outs)
                        outs.write('\n')
                        c += 1
                print('patients {}'.format(c))
    def get_patient_study():
        with tcia_client.get_patient_study() as response:
            if response.status == 200:
                data = response.read().decode('utf-8')
                get_patient_study = json.loads(data)
                c = 0
                with open('source/tcia/patient_study.json', 'w') as outs:
                    for p in get_patient_study:
                        json.dump(p, outs)
                        outs.write('\n')
                        c += 1
                print('patient_study {}'.format(c))

    def get_series():
        with tcia_client.get_series() as response:
            if response.status == 200:
                data = response.read().decode('utf-8')
                series = json.loads(data)
                c = 0
                with open('source/tcia/series.json', 'w') as outs:
                    for s in series:
                        json.dump(s, outs)
                        c += 1
                        outs.write('\n')
        print('series {}'.format(c))

    def get_sop_uids():
        series_instance_uids = []
        with open('source/tcia/series.json', 'r') as ins:
            for line in ins:
                series_instance_uids.append(json.loads(line)['SeriesInstanceUID'])

        with open('source/tcia/sop_uids.json', 'w') as outs:
            c = 0
            for series_instance_uid in series_instance_uids:
                with tcia_client.get_sop_uids(SeriesInstanceUID=series_instance_uid) as response:
                    if response.status == 200:
                        data = response.read().decode('utf-8')
                        sop_uids = json.loads(data)
                        s = {'SeriesInstanceUID': series_instance_uid, 'sop_instance_uids': [s['sop_instance_uid'] for s in sop_uids]}
                        json.dump(s, outs)
                        c += 1
                        outs.write('\n')
        print('sop_uids {}'.format(c))



    get_patient()
    get_patient_study()
    get_series()
    get_sop_uids()


if __name__ == '__main__':
    log_fmt = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    pathlib.Path("source/tcia").mkdir(parents=True, exist_ok=True)
    logging.basicConfig(level=logging.DEBUG, format=log_fmt)

    # not used in this stub but often useful for finding various files
    project_dir = Path(__file__).resolve().parents[2]

    # find .env automagically by walking up directories until it's found, then
    # load up the .env entries as environment variables
    load_dotenv(find_dotenv())

    transform()
