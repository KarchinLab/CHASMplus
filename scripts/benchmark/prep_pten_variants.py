"""
File: prep_pten_variants.py
Author: Collin Tokheim
Email: ctokheim@jhu.edu
Github: ctokheim
Description: Formats PTEN variants from saturation mutagenesis into a useable format
"""
import pandas as pd
import argparse


def parse_arguments():
    info = 'Formats PTEN variants from saturation mutagenesis into a useable format'
    parser = argparse.ArgumentParser(description=info)
    parser.add_argument('-i', '--input',
                        type=str, required=True,
                        help='saturation mutagenesis supp table for pten')
    parser.add_argument('-o', '--output',
                        type=str, required=True,
                        help='formatted output')
    args = parser.parse_args()
    return vars(args)


def main(opts):
    # read data
    df = pd.read_table(opts['input'])

    # add needed columns
    df['HGVSp_Short'] = 'p.' + df['Variant (one letter)']
    df['Hugo_Symbol'] = 'PTEN'

    # filter to high confidence missense
    is_missense = df['Type']=='missense'
    is_high_conf = df['High_conf']
    output_df = df[is_missense & is_high_conf]

    # save output
    output_df.to_csv(opts['output'], sep='\t', index=False)


if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)

