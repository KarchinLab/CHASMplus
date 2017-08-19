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
    coord_df = pd.read_table(opts['genomic_coord']).drop_duplicates()

    # keep only missense mutations
    df = df[df.aa.str.match('^[A-Z][0-9]+[A-Z]$')]

    # merge the info
    mut_info_df = coord_df['coordinates(gDNA/cDNA/protein)'].str.split('/', expand=True).rename(columns={ 0: 'gDNA', 1: 'cDNA', 2: 'protein'})
    mut_info_df['Chromosome'] = mut_info_df.gDNA.str.extract('^[chr]+([0-9XY]+)', expand=False)
    mut_info_df['Start_Position'] = mut_info_df.gDNA.str.extract(':g.([0-9]+)', expand=False)
    mut_info_df['End_Position'] = mut_info_df['Start_Position']
    mut_info_df['Reference_Allele'] = mut_info_df.gDNA.str.extract('[0-9]+([ACTG])>', expand=False)
    mut_info_df['Tumor_Seq_Allele2'] = mut_info_df.gDNA.str.extract('[0-9]+[ACTG]>([ACTG])$', expand=False)
    tmp = coord_df.input.str.split(':', expand=True).rename(columns={0: 'gene', 1: 'aa'})
    mut_info_df['gene'] = tmp['gene']
    mut_info_df['aa'] = tmp['aa']

    # merge in the information
    all_cols = ['Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele',
                'Tumor_Seq_Allele2', 'gene', 'aa']
    df = pd.merge(df, mut_info_df[all_cols], on=['gene', 'aa'], how='left')

    # Fix some of the columns
    df['HGVSp_Short'] = 'p.' + df['aa']
    df['Strand'] = '+'
    df['Variant_Classification'] = 'Missense_Mutation'
    df = df.rename(columns={'gene': 'Hugo_Symbol', 'consensus_high_call_DF0629': 'status'})

    # save merged results
    df.to_csv(opts['output'], sep='\t', index=False)


if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)

