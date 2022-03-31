#!/usr/bin/env python

import os
import logging
import argparse
from ftplib import FTP


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--output", default="source/pubmed/baseline")
    parser.add_argument("-N", default=1, type=int)
    args = parser.parse_args()

    ftp = FTP('ftp.ncbi.nlm.nih.gov')
    ftp.login()
    ftp.cwd('/pubmed/baseline/') 
    for i in ftp.nlst("/pubmed/baseline/*.xml.gz"):
        name = os.path.basename(i)
        localPath = os.path.join(args.output, name)
        if not os.path.exists( localPath ):
            print("download:", name)
            with open(localPath, 'wb') as fp:
                ftp.retrbinary('RETR %s' % (name), fp.write)
        else:
            print("cached:", name)
    ftp.quit()