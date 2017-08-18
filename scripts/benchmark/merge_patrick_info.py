"""
File: merge_patrick_info.py
Author: Collin Tokheim
Email: ctokheim@jhu.edu
Github: ctokheim
Description: Merge genomic coordinate information into Gordon mills assay
"""
import pandas as pd
import argparse


def parse_arguments():
    info = 'Merge genomic coordinate information into Gordon mills assay'
    parser = argparse.ArgumentParser(description=info)
    parser.add_argument('-i', '--input',
                        type=str, required=True,
                        help='Experiment file with functional classes')
    parser.add_argument('-g', '--genomic-coord',
                        type=str, required=True,
                        help='Experiment file with genomic coordinates')
    parser.add_argument('-o', '--output',
                        type=str, required=True,
                        help='Merge file')
    args = parser.parse_args()
    return vars(args)


def main(opts):
    # read data
    df = pd.read_table(opts['input'])
    coord_df = pd.read_table(opts['genomic_coord'])

    # merge the info
    mut_info_df = coord_df['coordinates(gDNA/cDNA/protein)'].str.split('/', expand=True).rename(columns={ 0: 'gDNA', 1: 'cDNA', 2: 'protein'})
    df['Chromosome'] = mut_info_df.gDNA.str.extract('^[chr]+([0-9]+)', expand=False)
    df['Start_Position'] = mut_info_df.gDNA.str.extract(':g.([0-9]+)', expand=False)
    df['Reference_Allele'] = mut_info_df.gDNA.str.extract('[0-9]+([ACTG])>', expand=False)
    df['Tumor_Seq_Allele2'] = mut_info_df.gDNA.str.extract('[0-9]+[ACTG]>([ACTG])$', expand=False)
    df['HGVSp_Short'] = 'p.' + df['aa']
    df = df.rename(columns={'gene': 'Hugo_Symbol', 'consensus_high_call_DF0629': 'status'})

    # save merged results
    df.to_csv(opts['output'], sep='\t', index=False)


if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)

