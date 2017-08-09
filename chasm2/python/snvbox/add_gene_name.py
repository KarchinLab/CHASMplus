"""
File: add_gene_name.py
Author: Collin Tokheim
Email: ctokheim@jhu.edu
Github: ctokheim
Description: Merges the gene name into snvbox output
"""
import pandas as pd
import csv
import argparse

# logging
import logging
logger = logging.getLogger(__name__)  # module logger

def parse_arguments():
    info = 'Add gene name to snvbox output'
    parser = argparse.ArgumentParser(description=info)
    parser.add_argument('-i', '--id2gene',
                        type=str, required=True,
                        help='file containing mapping of variant IDs to gene names')
    parser.add_argument('-s', '--snvbox',
                        type=str, required=True,
                        help='SNVbox features in tab delimited format')
    parser.add_argument('-f', '--features',
                        type=str, required=True,
                        help='Final output containing merged information')
    args = parser.parse_args()
    return vars(args)


def main(opts):
    logger.info('Adding gene names . . .')
    # read data
    id2gene_df = pd.read_table(opts['id2gene'], header=None, names=['UID', 'gene'])
    snvbox_df = pd.read_table(opts['snvbox'])

    # merge data
    merged_df = pd.merge(id2gene_df, snvbox_df, on='UID')

    # save results
    if opts['features']:
        merged_df.to_csv(opts['features'], sep='\t', index=False)
    logger.info('Done adding gene names.')
    return merged_df


if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)
