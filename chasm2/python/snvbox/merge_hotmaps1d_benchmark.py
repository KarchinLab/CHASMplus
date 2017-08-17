"""
File: merge_hotmaps1d_precompute.py
Author: Collin Tokheim
Email: ctokheim@jhu.edu
Github: ctokheim
Description: Merge the hotmaps1d features into snvbox features
"""
import pandas as pd
import argparse

# logging
import logging
logger = logging.getLogger(__name__)  # module logger


def parse_arguments():
    info = 'Merge the hotmaps1d features into snvbox features'
    parser = argparse.ArgumentParser(description=info)
    parser.add_argument('-hm', '--hotmaps',
                        type=str, required=True,
                        help='Hotmaps1d result file')
    parser.add_argument('-f', '--features',
                        type=str, required=True,
                        help='Features file from snvget')
    parser.add_argument('-o', '--output',
                        type=str, required=True,
                        help='Feature file with merged hotmaps result')
    args = parser.parse_args()
    return vars(args)


def main(opts):
    # merge the two data frames
    logger.info('merging hotmaps. . .')
    # read hotmaps
    if type(opts['hotmaps']) is str:
        hotmaps_df = pd.read_table(opts['hotmaps'])
    else:
        hotmaps_df = opts['hotmaps']

    # save class series
    logger.info('merge previous features ...')
    if type(opts['features']) is str:
        df = pd.read_table(opts['features'])
    else:
        df = opts['features']
    class_ser = df['class']
    del df['class']
    mycols = ['UID', 'gene',
              'p-value (0)', 'p-value (5)', 'p-value (10)']
    final_df = pd.merge(df, hotmaps_df.rename(columns={'index':'UID',
                                                      'Hugo_Symbol': 'gene'
                                                      })[mycols],
                        on='UID', how='left')
    final_df['class'] = class_ser
    if opts['output']:
        final_df.to_csv(opts['output'], sep='\t', index=False)

    return final_df


if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)
