"""
File: create_hg38_maf.py
Author: Collin Tokheim
Email: ctokheim@jhu.edu
Github: ctokheim
Description: Create a hg38 MAF from the liftover results
"""
import pandas as pd
import argparse


def parse_arguments():
    info = 'Create a hg38 MAF from the liftover results'
    parser = argparse.ArgumentParser(description=info)
    parser.add_argument('-i', '--input',
                        type=str, required=True,
                        help='Bed file from liftover')
    parser.add_argument('-m', '--maf',
                        type=str, required=True,
                        help='MAF file in hg19 coordinates')
    parser.add_argument('-o', '--output',
                        type=str, required=True,
                        help='MAF file with hg38 coordinates')
    args = parser.parse_args()
    return vars(args)


def main(opts):
    # read in bed file from liftover
    col_names = ['Chromosome_hg38', 'Start_Position_hg38',
                 'End_Position_hg38', 'index']
    hg38_bed = pd.read_table(opts['input'], header=None, names=col_names)
    hg38_bed['Start_Position_hg38'] = hg38_bed['Start_Position_hg38'] + 1
    hg38_bed['Chromosome_hg38'] = hg38_bed['Chromosome_hg38'].str[3:]

    # read in maf
    mut_df = pd.read_table(opts['maf'])
    mut_df['index'] = list(range(len(mut_df)))
    mut_df = mut_df.rename(columns={'Protein_Change': 'HGVSp_Short',
                                    'Tumor_Allele': 'Tumor_Seq_Allele2',
                                    'Gene': 'Hugo_Symbol'})
    column_order = mut_df.columns.tolist()

    # merge in bed information
    merged_df = pd.merge(mut_df, hg38_bed, on='index', how='left')

    # remove old coordinates and rename new ones
    del merged_df['Chromosome']
    del merged_df['Start_Position']
    del merged_df['End_Position']
    rename_dict = {
        'Chromosome_hg38': 'Chromosome',
        'Start_Position_hg38': 'Start_Position',
        'End_Position_hg38': 'End_Position',
    }
    merged_df = merged_df.rename(columns=rename_dict)

    # convert start/end positions back to int
    #merged_df = merged_df.dropna()
    merged_df['Chromosome'] = merged_df['Chromosome'].fillna('chrLiftoverFail')
    merged_df['Start_Position'] = merged_df['Start_Position'].fillna(-1)
    merged_df['End_Position'] = merged_df['End_Position'].fillna(-1)
    merged_df['Start_Position'] = merged_df['Start_Position'].astype(int)
    merged_df['End_Position'] = merged_df['End_Position'].astype(int)

    # save new maf file
    merged_df[column_order].to_csv(opts['output'], sep='\t', index=False)


if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)

