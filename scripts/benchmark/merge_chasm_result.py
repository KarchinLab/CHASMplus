"""
File: merge_chasm_result.py
Author: Collin Tokheim
Email: ctokheim@jhu.edu
Github: ctokheim
Description: Merge the CHASM result
"""
import os
import pandas as pd
import argparse


def parse_arguments():
    info = 'Merge the chasm results'
    parser = argparse.ArgumentParser(description=info)
    parser.add_argument('-i', '--input-dir',
                        type=str, required=True,
                        help='Directory containing files, including file prefix')
    parser.add_argument('-o', '--output',
                        type=str, required=True,
                        help='output file')
    args = parser.parse_args()
    return vars(args)


def main(opts):
    useful_cols = ['Hugo_Symbol',
                   #'Transcript_ID',
                   'HGVSp_Short',
                   'ID',
                   'class']
    # pattern = os.path.join(opts['input_dir'], 'driver_snvbox_inputdriver*.txt')
    df_list = []
    for i in range(10):
        path = opts['input_dir'] + '_chasm_input{0}driver{0}.output'.format(i)
        chasm_result = pd.read_table(path)
        path = opts['input_dir'] + '_split_driver{0}.txt'.format(i)
        mut_df = pd.read_table(path, usecols=useful_cols)
        #mut_df['MutationID'] = range(len(mut_df))
        merged = pd.merge(chasm_result, mut_df,
                          left_on='MutationID', right_on='ID',
                          how='left')
        df_list.append(merged)
    # get the passenger data
    path = opts['input_dir'] + '_chasm_input_passengerdriver0.output'
    chasm_result = pd.read_table(path)
    path = opts['input_dir'] + '_split_passenger.txt'.format(i)
    mut_df = pd.read_table(path, usecols=useful_cols)
    #mut_df['MutationID'] = range(len(mut_df))
    merged = pd.merge(chasm_result, mut_df,
                      left_on='MutationID', right_on='ID',
                      how='left')
    df_list.append(merged)

    # save results
    full_df = pd.concat(df_list)
    full_df.to_csv(opts['output'], sep='\t', index=False)


if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)


