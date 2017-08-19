"""
File: maf_to_annovar.py
Author: Collin Tokheim
Email: ctokheim@jhu.edu
Github: ctokheim
Description: Convert MAF file to input for ANNOVAR
"""
import pandas as pd
import argparse


def parse_arguments():
    info = 'Convert maf file to input for ANNOVAR'
    parser = argparse.ArgumentParser(description=info)
    parser.add_argument('-i', '--input',
                        type=str, required=True,
                        help='MAF file')
    parser.add_argument('-o', '--output',
                        type=str, required=True,
                        help='Path to annovar format file')
    args = parser.parse_args()
    return vars(args)


def main(opts):
    df = pd.read_table(opts['input'])
    out_cols = ['Chromosome', 'Start_Position', 'End_Position',
                'Reference_Allele', 'Tumor_Seq_Allele2']
    df[out_cols].to_csv(opts['output'], sep='\t', index=False, header=None)


if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)

