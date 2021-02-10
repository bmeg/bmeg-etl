import logging
import glob
import os
import ujson
import shutil

from bmeg.ioutils import reader
from bmeg.util.cli import default_argument_parser
from bmeg.util.logging import default_logging
from bmeg.emitter import new_emitter
from bmeg.stores import new_store
from bmeg import (Allele, SomaticCallset, SomaticCallset_Alleles_Allele)


IGNORE = "ignore"
TRANSLATE = "translate"
COPY = "copy"

opMapping = {
    "Allele.Vertex.json.gz" : IGNORE,
    "SomaticCallset.Vertex.json.gz" : COPY,
    "Aliquot_SomaticCallsets_SomaticCallset.Edge.json.gz" : COPY,
    "SomaticCallset_Aliquots_Aliquot.Edge.json.gz" : COPY,
    "SomaticCallset_Alleles_Allele.Edge.json.gz" : TRANSLATE,
    "Allele_SomaticCallsets_SomaticCallset.Edge.json.gz" : IGNORE
}

def init(
         store_path="outputs/allele/sqlite.db",
         annotated_alleles='outputs/allele/normalized.Allele.Vertex.json.gz'):
    logging.info('LOADING NORMALIZED ALLELES')
    logging.info("STORE PATH: %s", store_path)
    if os.path.isfile(store_path):
        os.remove(store_path)
    store = new_store('key-val', path=store_path, index=True)
    ac = 0
    with reader(annotated_alleles) as ins:
        for line in ins:
            vertex = ujson.loads(line)
            store.put(vertex['gid'], vertex['data'])
            ac += 1
            if ac % 100000 == 0:
                logging.info('alleles read: {}'.format(ac))

    logging.info('DONE')    
        
def transform(
              project_prefix,
              store_path="outputs/allele/sqlite.db",
              emitter_directory='allele',
              emitter_prefix='normalized',
              emitter_name='json'):

    logging.info("Scanning %s" % (project_prefix))
    filenames = []
    for filename in glob.iglob(project_prefix + "*"):
            filenames.append(filename)

    store = new_store('key-val', path=store_path, index=True)

    logging.info('FILES: {}'.format(filenames))
    ec = 0
    for f in filenames:
        op = IGNORE
        for k,v in opMapping.items():
            if f.endswith(k):
                op = v

        if op == TRANSLATE:
            logging.info('READING: {}'.format(f))
            
            input_directory = os.path.basename(os.path.dirname(f))
            input_prefix = os.path.basename(f).split(".")[0]
            emitter = new_emitter(name=emitter_name, directory=input_directory, prefix=input_prefix)
            with reader(f) as ins:
                for line in ins:
                    edge = ujson.loads(line)
                    data = store.get(edge['to'])
                    if data:
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
                    else:
                        logging.error("no allele found for: {}".format(edge))

                    if ec % 100000 == 0:
                        logging.debug('edges emitted: {}'.format(ec))
                emitter.close()
        elif op == COPY:
            input_directory = os.path.basename(os.path.dirname(f))
            outdir = os.path.join("outputs", input_directory)
            print("Copy %s to %s" % (f, outdir))
            shutil.copy(f, outdir)
        elif op == IGNORE:
            print("Ignoring file: %s" % (f))
        else:
            print("Unknown file: %s" % (f))
            
    logging.debug('total edges emitted: {}'.format(ec))



if __name__ == '__main__':  # pragma: no cover
    parser = default_argument_parser()
    parser.add_argument("--init", action="store_true", default=False)
    parser.add_argument('--project-prefix', type=str)
    options = parser.parse_args()
    default_logging(options.loglevel)
    if options.init:
        init()
    else:
        transform(project_prefix=options.project_prefix)
