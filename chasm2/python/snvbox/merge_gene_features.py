"""
File: merge_gene_features.py
Author: Collin Tokheim
Email: ctokheim@jhu.edu
Github: ctokheim
Description: Merge additional gene-level features into feature file
"""
import pandas as pd
import argparse

# logging
import logging
logger = logging.getLogger(__name__)  # module logger

def parse_arguments():
    info = 'Merge additional gene-level features into feature file'
    parser = argparse.ArgumentParser(description=info)
    parser.add_argument('-f', '--features',
                        type=str, required=True,
                        help='Feature file from snvbox (tab-delimited)')
    parser.add_argument('-d', '--dgd',
                        type=str, required=True,
                        help='File from duplicate gene database')
    parser.add_argument('-gl', '--gene-length',
                        type=str, required=True,
                        help='File with CDS length of gene')
    parser.add_argument('-l', '--lawrence',
                        type=str, required=True,
                        help='File containing HiC and replication timing from lawrence et al.')
    parser.add_argument('-o', '--output',
                        type=str, required=True,
                        help='Output feature file with merged features')
    args = parser.parse_args()
    return vars(args)


def main(opts):
    logger.info('merging gene features . . .')
    # read in data
    if opts['features'] is str:
        feat_df = pd.read_table(opts['features'])
    else:
        # if not file path, then it's a data frame
        feat_df = opts['features']
    dgd_df = pd.read_table(opts['dgd'])
    lawrence_df = pd.read_table(opts['lawrence'])
    gene_len_df = pd.read_table(opts['gene_length'])

    # modify columns for dgd data
    dgd_keep_cols = ['gene', 'num_dup_genes']
    rename_dict = {'Name': 'gene', 'NB_Genes': 'num_dup_genes'}
    # for some reason DGD has multiple rows for the same gene
    # so I have applied drop_duplicates
    dgd_df = dgd_df.rename(columns=rename_dict).drop_duplicates('gene')[dgd_keep_cols]

    # temporarily remove class column
    class_col = feat_df['class']  # save class column
    del feat_df['class']  # temporarily remove column

    # merge together data
    merged_df = pd.merge(feat_df, dgd_df, on='gene', how='left')
    lawrence_cols = ['gene', 'replication_time', 'HiC_compartment']
    merged_df = pd.merge(merged_df, lawrence_df[lawrence_cols], on='gene', how='left')
    tmp_cols = ['gene', 'gene length']
    merged_df = pd.merge(merged_df, gene_len_df[tmp_cols], on='gene', how='left')

    # perform appropriate fill of missing data
    merged_df['num_dup_genes'] = merged_df['num_dup_genes'].fillna(0)
    merged_df['replication_time'] = merged_df['replication_time'].fillna(merged_df['replication_time'].mean())
    merged_df['HiC_compartment'] = merged_df['HiC_compartment'].fillna(merged_df['HiC_compartment'].mean())
    merged_df['gene length'] = merged_df['gene length'].fillna(merged_df['gene length'].mean())

    # add back class column as last column
    merged_df['class'] = class_col

    # save output
    if opts['output']:
        merged_df.to_csv(opts['output'], sep='\t', index=False)
    logger.info('Done merging gene features.')
    return merged_df


if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)


