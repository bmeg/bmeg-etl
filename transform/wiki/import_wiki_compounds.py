import json
import requests
import click
from pprint import pprint

# gen3 interface ----------------------------------
from gen3.auth import Gen3Auth


def make_auth(endpoint, refresh_file):
    """Wraps Gen3Auth."""
    return Gen3Auth(endpoint, refresh_file=refresh_file)


def validate(auth, endpoint, path):
    """Tests endpoint/path for 200 response"""
    r = requests.get(f"{endpoint}/{path}", auth=auth)
    assert r.status_code == 200, f"Should return 200 {r}"


printouts = '\n'.join("""?Term=term
?GID=gid
?SelectedOntology.TermId=selected_ontology_term_id
?SelectedOntology.Term=selected_ontology_term
?Ontology.TermId=ontology_term_id
?Ontology.Term=ontology_term
?CompoundSynonym.Term=synonym_terms
""".split())

limit = 100

offset = 0

# see https://github.com/SemanticMediaWiki/SemanticMediaWiki/issues/3168
params = {
    'title': 'Special:Ask',
    'q': '[[Category:Compound]]',
    'po': printouts,
    'p[format]': 'json',
    'p[limit]': limit,
    'p[offset]': offset,
    'request_type': 'raw',
    'format': 'json'
}


@click.command()
@click.option('--endpoint', default='https://gen3.compbio.ohsu.edu', help='gen3 endpoint [https://gen3.compbio.ohsu.edu]')
@click.option('--path', default='bmeg-wiki/index.php', help='wiki endpoint [bmeg-wiki/index.php]')
@click.option('--output', default='outputs/wiki/CompoundSynonyms.json', help='output file path')
@click.option('--refresh_file', default='credentials.json', help='gen3 credentials file')
def transform(endpoint, path, output, refresh_file):
    print(endpoint, path, refresh_file, output)
    auth = make_auth(endpoint, refresh_file)
    working = True
    with open(output, 'w') as out_fh:
        while working:
            url = f"{endpoint}/{path}"
            r = requests.post(url, data=params, auth=auth)
            try:
                assert r.status_code == 200, f"Should return 200 {r} {url} {r.text}"
                response = r.json()
            except Exception as e:
                print(f"ERROR:\nparams:{params}\nresponse was: {r.text} exception: {str(e)}")
                raise
            rows = response['rows']
            if rows < limit:
                working = False
            params['p[offset]'] += limit
            results = response['results']
            for compound_key, compound in results.items():
                po = compound['printouts']
                lookup = {
                    'compound_gid': po['gid'][0],
                    'selected_ontology': {'term': po['selected_ontology_term'][0], 'term_id': po['selected_ontology_term_id'][0]},
                    'ontologies': [{'term': o[0], 'term_id': o[1]} for o in zip(po['ontology_term'], po['ontology_term_id'])],
                    'synonyms': [s for s in po['term'] + po['synonym_terms']]
                }
                pprint(compound_key)
                pprint(po)
                pprint(lookup)
                exit()
                json.dump(lookup, out_fh, separators=(',', ':'))
                out_fh.write('\n')
    print(f'Done {output}')


if __name__ == '__main__':
    transform()
