"""
File: compute_feature_importance.py
Author: Collin Tokheim
Email: ctokheim@jhu.edu
Github: ctokheim
Description: Computes feature importance based on permutation
"""
import numpy as np
import pandas as pd
import argparse
import glob
import os


def cummin(x):
    """A python implementation of the cummin function in R"""
    for i in range(1, len(x)):
        if x[i-1] < x[i]:
            x[i] = x[i-1]
    return x


def bh_fdr(pval):
    """A python implementation of the Benjamani-Hochberg FDR method.
    This code should always give precisely the same answer as using
    p.adjust(pval, method="BH") in R.
    Parameters
    ----------
    pval : list or array
        list/array of p-values
    Returns
    -------
    pval_adj : np.array
        adjusted p-values according the benjamani-hochberg method
    """
    pval_array = np.array(pval)
    sorted_order = np.argsort(pval_array)
    original_order = np.argsort(sorted_order)
    pval_array = pval_array[sorted_order]

    # calculate the needed alpha
    n = float(len(pval))
    pval_adj = np.zeros(int(n))
    i = np.arange(1, int(n)+1, dtype=float)[::-1]  # largest to smallest
    pval_adj = np.minimum(1, cummin(n/i * pval_array[::-1]))[::-1]
    return pval_adj[original_order]


def parse_arguments():
    info = 'Computes feature importance based on permutation'
    parser = argparse.ArgumentParser(description=info)
    parser.add_argument('-i', '--input',
                        type=str, required=True,
                        help='Feature importance of observed training data')
    parser.add_argument('-p', '--permutation',
                        type=str, required=True,
                        help='Directory containing feature importance on permuted results')
    parser.add_argument('-o', '--output',
                        type=str, required=True,
                        help='output result file')
    args = parser.parse_args()
    return vars(args)


def main(opts):
    # read in actual feature importance
    df = pd.read_table(opts['input'])

    # create data frame of permuted feature importances
    df_list = []
    for f in glob.glob(os.path.join(opts['permutation'], 'variable_importance[0-9]*.txt')):
        tmp_df = pd.read_table(f)
        df_list.append(tmp_df)
    permuted_df = pd.concat(df_list, axis=1)

    # calculate p-value
    num_perm = len(permuted_df.columns)
    pval = permuted_df.apply(lambda x: x>=df['MeanDecreaseGini'], axis=0).sum(axis=1) / float(num_perm)
    output_df = pval.reset_index(name='pvalue').rename(columns={'index': 'feature'})

    # calculate fdr
    output_df['qvalue'] = bh_fdr(output_df['pvalue'])

    # add mean decrease gini
    output_df['MeanDecreaseGini'] = df['MeanDecreaseGini'].tolist()

    # calculate z-score
    mean = permuted_df.mean(axis=1)
    std = permuted_df.std(axis=1)
    zscore = ((output_df.set_index('feature')['MeanDecreaseGini'] - mean) / std)
    output_df['zscore'] = zscore.tolist()

    # save results
    output_df.to_csv(opts['output'], sep='\t', index=False)


if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)

