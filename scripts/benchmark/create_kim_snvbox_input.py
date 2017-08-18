"""
File: create_kim_snvbox_input.py
Author: Collin Tokheim
Email: ctokheim@jhu.edu
Github: ctokheim
Description: Create transcript input format for snvbox for kim et al data
"""
import pandas as pd
import argparse
import os
import glob


def parse_arguments():
    info = 'Create transcript input format for snvbox for kim et al data'
    parser = argparse.ArgumentParser(description=info)
    parser.add_argument('-i', '--input',
                        type=str, required=True,
                        help='Kim et al mutation data')
    parser.add_argument('-p', '--preferred-tx',
                        type=str, required=True,
                        help='File containing the preferred transcript for each gene')
    parser.add_argument('-o', '--output',
                        type=str, required=True,
                        help='Snvbox output file (or dir if cross-val specified)')
    args = parser.parse_args()
    return vars(args)


def main(opts):
    # read kim et al mutations
    df = pd.read_table(opts['input'])
    #df = df.Mutation.str.split('_', expand=True).rename(columns={0: 'gene', 1: 'mutation'})
    df['mutation'] = df.mutation.str[2:]

    # read preferred snvbox transcripts
    tx_df = pd.read_table(opts['preferred_tx'])

    # merge the two
    merged_df = pd.merge(df, tx_df, on='gene')

    # add a UID column
    merged_df['UID'] = range(len(merged_df))

    # save results
    mycols = ['UID', 'transcript', 'mutation']
    merged_df[mycols].to_csv(opts['output'], sep=' ', index=False, header=None)


if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)
