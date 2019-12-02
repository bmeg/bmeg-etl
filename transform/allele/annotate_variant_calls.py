import logging
import glob
import ujson
from bmeg.ioutils import reader

from bmeg.util.cli import default_argument_parser
from bmeg.util.logging import default_logging
from bmeg.emitter import new_emitter
from bmeg import (Allele, SomaticCallset, SomaticCallset_Alleles_Allele)


def transform(output_dir='outputs',
              edge_filename_pattern='**/*SomaticCallset_Alleles_Allele.Edge.json.gz',
              annotated_alleles='outputs/allele/normalized.Allele.Vertex.json.gz',
              emitter_directory='allele',
              emitter_prefix='normalized',
              emitter_name='json'):

    emitter = new_emitter(name=emitter_name, directory=emitter_directory, prefix=emitter_prefix)

    path = '{}/{}'.format(output_dir, edge_filename_pattern)
    filenames = []
    for filename in glob.iglob(path, recursive=True):
        if 'outputs/allele/' in filename:
            continue
        else:
            filenames.append(filename)

    logging.info('LOADING NORMALIZED ALLELES')
    alleles = {}
    with reader(annotated_alleles) as ins:
        for line in ins:
            vertex = ujson.loads(line)
            alleles[vertex['gid']] = vertex['data']
    logging.info('DONE')

    logging.info('FILES: {}'.format(filenames))
    ec = 0
    for f in filenames:
        logging.info('READING: {}'.format(f))
        with reader(f) as ins:
            for line in ins:
                try:
                    edge = ujson.loads(line)
                    if edge['to'] in alleles:
                        data = alleles[edge['to']]
                        edge['data']['ensembl_gene'] = data.get('ensembl_gene')
                        edge['data']['ensembl_transcript'] = data.get('ensembl_transcript')
                        edge['data']['ensembl_protein'] = data.get('ensembl_protein')
                        allele_gid = Allele.make_gid(data['genome'],
                                                     data['chromosome'],
                                                     data['start'],
                                                     data['reference_bases'],
                                                     data['alternate_bases'])
                        assert allele_gid == edge['to']
                        emitter.emit_edge(
                            SomaticCallset_Alleles_Allele(
                                from_gid=SomaticCallset.make_gid(*edge['from'].replace('SomaticCallset:', '').split(':')),
                                to_gid=allele_gid,
                                data=edge['data']
                            ),
                            emit_backref=True
                        )
                        ec += 1
                except Exception as e:
                    logging.error(str(e))

                if ec % 10000 == 0:
                    logging.debug('edges emitted: {}'.format(ec))

    logging.debug('total edges emitted: {}'.format(ec))
    emitter.close()


if __name__ == '__main__':  # pragma: no cover
    parser = default_argument_parser()
    options = parser.parse_args()
    default_logging(options.loglevel)
    transform()
