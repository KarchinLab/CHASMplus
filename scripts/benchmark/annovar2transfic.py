"""
File: annovar2transfic.py
Author: Collin Tokheim
Email: ctokheim@jhu.edu
Github: ctokheim
Description: Convert annovar format to TransFIC format
"""
import pandas as pd
import argparse
import numpy as np


def parse_arguments():
    info = 'Convert annovar format to TransFIC format'
    parser = argparse.ArgumentParser(description=info)
    parser.add_argument('-i', '--input',
                        type=str, required=True,
                        help='Annotation file from annovar using ensemble genes')
    parser.add_argument('-o', '--output',
                        type=str, required=True,
                        help='TransFIC file')
    args = parser.parse_args()
    return vars(args)


def main(opts):
    # read data
    df = pd.read_table(opts['input'])
    df['mygene'] = df['Gene.ensGene'].str.extract('([A-Z0-9]+);*')
    df['mutation id'] = range(len(df))
    mycols = ['mutation id', 'mygene', 'SIFT_score', 'Polyphen2_HDIV_score', 'MutationAssessor_score']

    # fill in NAs
    is_na = df['SIFT_score'].astype(str)=='.'
    df.loc[is_na, 'SIFT_score'] = np.nan
    is_na = df['Polyphen2_HDIV_score'].astype(str)=='.'
    df.loc[is_na, 'Polyphen2_HDIV_score'] = np.nan
    is_na = df['MutationAssessor_score'].astype(str)=='.'
    df.loc[is_na, 'MutationAssessor_score'] = np.nan

    # drop columns with all na
    df = df[df[mycols].isnull().sum(axis=1)==0]

    # save results
    df[mycols].to_csv(opts['output'], sep='\t', index=False, header=None)


if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)


