#!/usr/bin/env python3

import argparse
import datetime
import os
import pandas
import requests
import sys

from urllib.parse import urlparse


def update_manifest(manifest, fname, url, description):
    print("Updating ", manifest, "...", file=sys.stderr, sep="")
    manifest_df = pandas.read_csv(manifest, sep="\t")
    now = datetime.datetime.now()
    new_row = pandas.DataFrame(data={"Name": [fname],
                                     "URL": [url],
                                     "Date": [now.strftime("%Y-%m-%d")],
                                     "Description": [description]})
    manifest_df = manifest_df.append(new_row)
    manifest_df.to_csv(manifest, sep="\t", index=False)


def download(url, fname):
    print("Downloading file: ", url, "...", file=sys.stderr, sep="")

    outdir = os.path.dirname(fname)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    elif not os.path.isdir(outdir):
        raise ValueError("%s exists and is not a directory" % outdir)

    response = requests.get(url)
    response.raise_for_status()

    with open(fname, 'wb') as handle:
        for block in response.iter_content(1024):
            handle.write(block)

    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('url', type=str,
                        help='Url of file to download')
    parser.add_argument('--output-dir', type=str, default=".",
                        help='Output directory')
    parser.add_argument('--name', type=str, required=False,
                        help='Name to use for downloaded file')
    parser.add_argument('--description', type=str, required=True,
                        help='Description to use for downloaded file')

    args = parser.parse_args()
    if not args.name:
        a = urlparse(args.url)
        args.name = os.path.basename(a.path)

    outdir = os.path.abspath(".")
    fname = os.path.join(outdir, args.name)

    download(args.url, fname)

    manifest = os.path.join(outdir, "manifest.tsv")
    update_manifest(manifest, args.name, args.url, args.description)
