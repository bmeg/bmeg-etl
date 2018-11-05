from bmeg.enrichers.transcript_enricher import normalize
from bmeg.emitter import JSONEmitter
import concurrent.futures

""" read a file of transcript_ids, call ensembl and emit vertex"""


def transform(output_dir,
              prefix,
              missing_transcript_ids_filename):
    emitter = JSONEmitter(output_dir, prefix)

    # threaded worker
    def worker(transcript_id):
        transcript = normalize(transcript_id)
        if transcript:
            emitter.emit_vertex(transcript)
    # main
    with open(missing_transcript_ids_filename) as ins:
        with concurrent.futures.ThreadPoolExecutor(8) as executor:
            for line in ins:
                transcript_id = line.strip()
                executor.submit(worker, transcript_id)
    emitter.close()


def main():  # pragma: no cover
    transform(output_dir='ensembl',
              prefix='missing',
              missing_transcript_ids_filename='source/ensembl/missing_transcript_ids.txt')


if __name__ == '__main__':  # pragma: no cover
    main()
