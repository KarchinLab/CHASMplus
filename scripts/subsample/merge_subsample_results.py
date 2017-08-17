"""
File: merge_subsample_results.py
Author: Collin Tokheim
Email: ctokheim@jhu.edu
Github: ctokheim
Description: Merge the subsample results into a summary
"""
import pandas as pd
import numpy as np
import os
import argparse


def parse_arguments():
    info = 'Merge the subsample results into a summary'
    parser = argparse.ArgumentParser(description=info)
    parser.add_argument('-i', '--input',
                        type=str, required=True,
                        help='Directory containing subsample results')
    parser.add_argument('-m', '--mutations',
                        type=str, required=True,
                        help='MAF file containing all mutations for all samples')
    parser.add_argument('-o', '--output',
                        type=str, required=True,
                        help='Output file')
    args = parser.parse_args()
    return vars(args)


def main(opts):
    # read data
    mypath = os.path.join(opts['input'], 'sample_rate.txt')
    samp_rate_df = pd.read_table(mypath)
    mut_df = pd.read_table(opts['mutations'])

    # loop through each iteration
    num_mut_list = []
    avg_mut_list = []
    for i in range(30):
        mypath = os.path.join(opts['input'],
                              'sample{0}/chasm2_result_pretrained_final.txt'.format(i))
        if os.path.exists(mypath):
            # read results and mutations
            df = pd.read_table(mypath)
            mypath = os.path.join(opts['input'], 'sample{0}/mutations.maf'.format(i))
            tmp_mut_df = pd.read_table(mypath)
            tmp_mut_df['UID'] = range(len(tmp_mut_df))

            # add the number of significant mutations
            signif = df[df['CHASM2_genome_qval']<.01].copy()
            num_mut_list.append(len(signif))

            # average mutations per sample
            mycols = ['UID', 'Hugo_Symbol', 'Transcript_ID', 'HGVSp_Short']
            merged = pd.merge(signif, tmp_mut_df[mycols], on='UID', how='left')
            merged['flag'] = 1
            merged = pd.merge(mut_df, merged,
                              on=['Hugo_Symbol', 'Transcript_ID', 'HGVSp_Short'],
                              how='left')
            avg = merged['flag'].sum() / merged.Tumor_Sample_Barcode.nunique()
            avg_mut_list.append(avg)
        else:
            num_mut_list.append(np.nan)
            avg_mut_list.append(np.nan)

    # add a new column with num mutations significant
    samp_rate_df['number signif mutations'] = num_mut_list
    samp_rate_df['avg mutations per sample'] = avg_mut_list

    # save results
    samp_rate_df.to_csv(opts['output'], sep='\t', index=False)


if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)

