"""
File: bed_to_snvbox.py
Author: Collin Tokheim
Email: ctokheim@jhu.edu
Github: ctokheim
Description:
"""
import pandas as pd
import argparse


def parse_arguments():
    info = 'Convert a liftovered bed file back to snvbox format'
    parser = argparse.ArgumentParser(description=info)
    parser.add_argument('-i', '--input',
                        type=str, required=True,
                        help='snvbox file before liftove')
    parser.add_argument('-b', '--bed',
                        type=str, required=True,
                        help='BED file')
    parser.add_argument('-o', '--output',
                        type=str, required=True,
                        help='snvbox file after liftover to save')
    args = parser.parse_args()
    return vars(args)


def main(opts):
    # read bed
    bed_names = ['chrom_new', 'start', 'pos_new', 'UID']
    bed = pd.read_table(opts['bed'], header=None, names=bed_names)

    # read snvbox input format
    snvbox_names = ['UID', 'chrom', 'pos', 'strand', 'ref', 'alt']
    snvbox = pd.read_table(opts['input'], header=None, names=snvbox_names)

    # merge
    merged = pd.merge(snvbox, bed, on='UID', how='left')

    # save file
    out_cols = ['UID', 'chrom_new', 'pos_new', 'strand', 'ref', 'alt']
    merged[out_cols].to_csv(opts['output'], sep='\t', index=False, header=None)


if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)


