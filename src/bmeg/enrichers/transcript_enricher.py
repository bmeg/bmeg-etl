from bmeg import Transcript, Project
import logging

from bmeg.requests import Client
requests = Client('transcript_enricher')


def normalize(ensembl_transcript_id):
    """ call ensembl and retrieve Transcript """
    url = 'https://rest.ensembl.org/lookup/id/{}?content-type=application/json;expand=1'.format(ensembl_transcript_id)
    r = requests.get(url, timeout=20)
    response = r.json()
    # {
    #   "db_type":"core",
    #   "start":55287849,
    #   "Parent":"ENSG00000115355",
    #   "end":55296455,
    #   "object_type":"Transcript",
    #   "logic_name":"havana",
    #   "assembly_name":"GRCh38",
    #   "species":"homo_sapiens",
    #   "biotype":"protein_coding",
    #   "seq_region_name":"2",
    #   "id":"ENST00000647547",
    #   "display_name":"CCDC88A-268",
    #   "source":"havana",
    #   "strand":-1,
    #   "is_canonical":0,
    #   "version":1
    # }
    if 'strand' not in response:
        logging.warning('no response for {}'.format(ensembl_transcript_id))
        logging.warning(response)
        return None
    strand = '+'
    if response['strand'] > 0:
        strand = '-'
    return Transcript(
        submitter_id=Transcript.make_gid(response['id']),
        transcript_id=response['id'],
        chromosome=response['seq_region_name'],
        start=response['start'],
        end=response['end'],
        strand=strand,
        biotype=response['biotype'],
        genome=response['assembly_name'],
        project_id=Project.make_gid("Reference")
    )
