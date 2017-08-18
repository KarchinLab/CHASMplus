"""
File: snvbox_to_bed.py
Author: Collin Tokheim
Email: ctokheim@jhu.edu
Github: ctokheim
Description:
"""
import pandas as pd
import argparse


def parse_arguments():
    info = 'Convert snvbox input format to BED'
    parser = argparse.ArgumentParser(description=info)
    parser.add_argument('-i', '--input',
                        type=str, required=True,
                        help='snvbox genomic file')
    parser.add_argument('-o', '--output',
                        type=str, required=True,
                        help='BED file')
    args = parser.parse_args()
    return vars(args)


def main(opts):
    # read data
    cols = ['UID', 'chrom', 'pos', 'strand', 'ref', 'alt']
    df = pd.read_table(opts['input'], header=None, names=cols)

    # zero based
    df['start0'] = df['pos'] - 1

    # save file
    mycols = ['chrom', 'start0', 'pos', 'UID']
    df[mycols].to_csv(opts['output'], sep='\t', index=False, header=None)


if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)

