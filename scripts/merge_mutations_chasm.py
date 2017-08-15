"""
File: merge_mutations.py
Author: Collin Tokheim
Email: ctokheim@jhu.edu
Github: ctokheim
Description: Merge mutations with CHASM result
"""
import pandas as pd
import argparse


def parse_arguments():
    info = 'Merge mutations with CHASM result'
    parser = argparse.ArgumentParser(description=info)
    parser.add_argument('-m', '--mutations',
                        type=str, required=True,
                        help='MAF file')
    parser.add_argument('-c', '--chasm',
                        type=str, required=True,
                        help='CHASM result file')
    parser.add_argument('-o', '--output',
                        type=str, required=True,
                        help='')
    args = parser.parse_args()
    return vars(args)


def main(opts):
    # read data
    mut_df = pd.read_table(opts['mutations'])
    mut_df['UID'] = range(len(mut_df))
    useful_cols = ['ID', 'UID', 'driver score', 'CHASM2', 'CHASM2_genome',
                   'CHASM2_pval', 'CHASM2_qval', 'CHASM2_genome_pval',
                   'CHASM2_genome_qval']
    chasm_df = pd.read_table(opts['chasm'], usecols=useful_cols)

    # merge data
    cols_of_interest = ['UID', 'Hugo_Symbol', 'Transcript_ID', 'HGVSp_Short']
    chasm_df = pd.merge(chasm_df, mut_df[cols_of_interest],
                        on='UID', how='left')
    cols_of_interest = ['Hugo_Symbol', 'Transcript_ID', 'HGVSp_Short',
                        'driver score', 'CHASM2', 'CHASM2_genome',
                        'CHASM2_pval', 'CHASM2_qval', 'CHASM2_genome_pval',
                        'CHASM2_genome_qval']
    mut_df = pd.merge(mut_df, chasm_df[cols_of_interest],
                      on=['Hugo_Symbol', 'Transcript_ID', 'HGVSp_Short'],
                      how='left')

    # save results
    mut_df.to_csv(opts['output'], sep='\t', index=False)


if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)

