"""
File: generate_chasm2_null_dist.py
Author: Collin Tokheim
Email: ctokheim@jhu.edu
Github: ctokheim
Description: Generate the empirical null distribution for CHASM2 scores
"""
import pandas as pd
import argparse
import os
import glob


def parse_arguments():
    info = 'Generate the empirical null distribution for CHASM2 scores'
    parser = argparse.ArgumentParser(description=info)
    parser.add_argument('-s', '--simulation-dir',
                        type=str, required=True,
                        help='Base simulation directory')
    parser.add_argument('-o', '--output',
                        type=str, required=True,
                        help='Output null distribution file')
    args = parser.parse_args()
    return vars(args)


def main(opts):
    sim_dir = opts['simulation_dir']

    # read in the 20/20+ results
    useful_cols = ['gene', 'driver score']
    tt_dict = {}
    for i in range(1, 11):
        path = os.path.join(sim_dir, '2020plus/sim{0}/results/r_random_forest_prediction.txt'.format(i))
        tt_dict[i] = pd.read_table(path, usecols=useful_cols)

    # go through the CHASM2 results
    useful_cols = ['gene', 'ID', 'UID', 'driver']
    df_list = []
    for i in range(1, 11):
        # read data
        path = os.path.join(sim_dir, 'output/chasm_result{0}.txt'.format(i))
        tmp_df = pd.read_table(path, usecols=useful_cols)

        # merge in 20/20+
        tmp_merged_df = pd.merge(tmp_df, tt_dict[i], on='gene', how='left')

        # create relevant scores
        tmp_merged_df['CHASM2'] = tmp_merged_df['driver']
        tmp_merged_df['CHASM2_genome'] = tmp_merged_df['driver'] * tmp_merged_df['driver score']

        # append to list
        df_list.append(tmp_merged_df[['CHASM2', 'CHASM2_genome']])

    # calculate empirical null
    cat_df = pd.concat(df_list)
    chasm2_gen_cts = cat_df['CHASM2_genome'].value_counts().sort_index(ascending=False)
    chasm2_mut_cts = cat_df['CHASM2'].value_counts().sort_index(ascending=False)
    full_index = chasm2_gen_cts.index.union(chasm2_mut_cts.index)
    output_df = pd.DataFrame({'CHASM2_genome_cts': chasm2_gen_cts,
                              'CHASM2_cts': chasm2_mut_cts},
                             index=full_index)
    output_df = output_df.sort_index(ascending=False)
    cum_sum_gen = output_df['CHASM2_genome_cts'].cumsum()
    cum_sum_mut = output_df['CHASM2_cts'].cumsum()
    output_df['CHASM2_genome'] = cum_sum_gen / output_df['CHASM2_genome_cts'].sum()
    output_df['CHASM2'] = cum_sum_mut / output_df['CHASM2_cts'].sum()

    # save results
    output_df.index.name = 'score'
    output_df = output_df.reset_index()
    out_cols = ['score', 'CHASM2', 'CHASM2_genome']
    output_df[out_cols].to_csv(opts['output'], sep='\t', index=False)


if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)


