"""
File: snvbox2candra.py
Author: Collin Tokheim
Email: ctokheim@jhu.edu
Github: ctokheim
Description: Convert snvbox genomic input into candra format
"""
import pandas as pd
import argparse


def parse_arguments():
    info = 'Convert snvbox genomic input into candra format'
    parser = argparse.ArgumentParser(description=info)
    parser.add_argument('-i', '--input',
                        type=str, required=True,
                        help='SNVbox genomic format')
    parser.add_argument('-o', '--output',
                        type=str, required=True,
                        help='Candra format')
    args = parser.parse_args()
    return vars(args)


def main(opts):
    # read data
    cols = ['mutation_id', 'chrom', 'start', 'strand', 'ref', 'alt']
    df = pd.read_table(opts['input'], header=None, names=cols)

    # fix chromosome
    df['chrom'] = df['chrom'].str[3:]

    # save results
    out_cols = ['chrom', 'start', 'ref', 'alt', 'strand']
    df[out_cols].to_csv(opts['output'], sep='\t', index=False, header=None)


if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)


