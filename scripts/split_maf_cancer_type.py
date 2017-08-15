"""
File: split_maf_cancer_type.py
Author: Collin Tokheim
Email: ctokheim@jhu.edu
Github: ctokheim
Description: Split the MAF file by cancer type
"""
import pandas as pd
import argparse
import os


def parse_arguments():
    info = 'Split the MAF file by cancer type'
    parser = argparse.ArgumentParser(description=info)
    parser.add_argument('-i', '--input',
                        type=str, required=True,
                        help='Full MAF file')
    parser.add_argument('-o', '--output',
                        type=str, required=True,
                        help='Output directory to save splits')
    args = parser.parse_args()
    return vars(args)


def main(opts):
    # read data
    df = pd.read_table(opts['input'])

    # save each cancer type
    cancer_types = df['CODE'].unique()
    for c in cancer_types:
        ctype_dir = os.path.join(opts['output'], c)
        if not os.path.exists(ctype_dir): os.makedirs(ctype_dir)
        tmp_path = os.path.join(ctype_dir, 'mutations.maf')
        df[df['CODE']==c].to_csv(tmp_path, sep='\t', index=False)


if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)

