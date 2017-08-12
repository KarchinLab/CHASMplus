"""
File: combined_score.py
Author: Collin Tokheim
Email: ctokheim@jhu.edu
Github: ctokheim
Description: Computes the combined CHASM2/20/20+ score and calculates p-values
"""
import pandas as pd
import chasm2.python.utils as utils
import argparse


def parse_arguments():
    info = 'Computes the combined score of 20/20+ and CHASM2, along with calculating p-values'
    parser = argparse.ArgumentParser(description=info)
    parser.add_argument('-nd', '--null-distribution',
                        type=str, required=True,
                        help='Null distribution file')
    parser.add_argument('-c', '--chasm2',
                        type=str, required=True,
                        help='Score results from CHASM2')
    parser.add_argument('-t', '--twentyTwentyPlus',
                        type=str, required=True,
                        help='Score results from 20/20+')
    parser.add_argument('-o', '--output',
                        type=str, required=True,
                        help='Final output from CHASM2')
    args = parser.parse_args()
    return vars(args)


def main(opts):
    # read results
    ttp_df = pd.read_table(opts['twentyTwentyPlus'])
    chasm_df = pd.read_table(opts['chasm2'])

    # merge of gene scores
    df = pd.merge(chasm_df, ttp_df[['gene', 'driver score', 'oncogene score', 'tsg score']],
                  on='gene', how='left')

    df['CHASM2_genome'] = df['driver score']*df['CHASM2']

    # read the empirical null distribution file
    null_dist = pd.read_table(opts['null_distribution'], index_col=0)

    # calc the p-values
    pvals_genome = utils.compute_p_value(df['CHASM2_genome'].dropna(), null_dist['CHASM2_genome'].dropna())
    pvals_regular = utils.compute_p_value(df['CHASM2'].dropna(), null_dist['CHASM2'].dropna())

    ## p-value/FDR for exome-wide
    # fill in the p-values and fdr
    df['CHASM2_genome_pval'] = np.nan
    df.loc[~df['CHASM2_genome'].isnull(), 'CHASM2_genome_pval'] = pvals_genome.tolist()
    df['CHASM2_genome_qval'] = np.nan
    df.loc[~df['CHASM2_genome'].isnull(), 'CHASM2_genome_qval'] = utils.bh_fdr(pvals_genome)

    ## p-value/FDR for regular CHASM2
    # fill in the p-values and fdr
    df['CHASM2_pval'] = np.nan
    df.loc[~df['CHASM2'].isnull(), 'CHASM2_pval'] = pvals_regular.tolist()
    df['CHASM2_qval'] = np.nan
    df.loc[~df['CHASM2'].isnull(), 'CHASM2_qval'] = utils.bh_fdr(pvals_regular)

    # save results
    df.to_csv(opts['output'], sep='\t', index=False)


if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)

