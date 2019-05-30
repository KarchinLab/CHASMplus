"""
File: maf_to_bed.py
Author: Collin Tokheim
Email: ctokheim@jhu.edu
Github: ctokheim
Description: Convert MAF file to BED file for use by liftOver
"""
import pandas as pd
import argparse


def parse_arguments():
    info = 'Convert MAF file to BED file for use by liftOver'
    parser = argparse.ArgumentParser(description=info)
    parser.add_argument('-i', '--input',
                        type=str, required=True,
                        help='MAF file')
    parser.add_argument('-o', '--output',
                        type=str, required=True,
                        help='Bed file')
    args = parser.parse_args()
    return vars(args)


def main(opts):
    # read maf
    useful_cols = ['Chromosome', 'Start_Position', 'End_Position']
    mut_df = pd.read_table(opts['input'], usecols=useful_cols)
    mut_df['name'] = list(range(len(mut_df)))

    # fix positions and chromosomes
    mut_df['Start_Position'] = mut_df['Start_Position'] - 1
    mut_df['Chromosome'] = mut_df['Chromosome'].astype(str)
    has_chr = mut_df['Chromosome'].str.startswith('chr')
    mut_df.loc[has_chr, 'Chromosome'] = mut_df.loc[has_chr, 'Chromosome'].str[3:]
    mut_df['Chromosome'] = 'chr' + mut_df['Chromosome'].astype(str)

    # save to BED file
    out_cols = ['Chromosome', 'Start_Position', 'End_Position', 'name']
    mut_df['End_Position'] = mut_df['End_Position'].astype(int)
    mut_df[out_cols].to_csv(opts['output'], sep='\t', index=False, header=None)


if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)

