"""
File: create_berger_snvbox_input.py
Author: Collin Tokheim
Email: ctokheim@jhu.edu
Github: ctokheim
Description: Create transcript input format for snvbox for IARC TP53 data
"""
import pandas as pd
import argparse
import os
import glob


def parse_arguments():
    info = 'Create transcript input format for snvbox for kim et al data'
    parser = argparse.ArgumentParser(description=info)
    parser.add_argument('-i', '--input',
                        type=str, required=True,
                        help='mutation data from berger et al')
    parser.add_argument('-p', '--preferred-tx',
                        type=str, required=True,
                        help='File containing the preferred transcript for each gene')
    parser.add_argument('-c', '--cross-val-dir',
                        type=str, default=None,
                        help='Cross-validation directory needed for original CHASM')
    parser.add_argument('-pf', '--prefix',
                        type=str, default=None,
                        help='File prefix for output')
    parser.add_argument('-o', '--output',
                        type=str, required=True,
                        help='Snvbox input format file')
    args = parser.parse_args()
    return vars(args)


def read_cross_val(mydir):
    # read training folds
    pattern = os.path.join(mydir, 'training*.txt')
    training_genes = dict()
    for f in glob.glob(pattern):
        fold = int(f[-5])
        with open(f) as handle:
            tmp_glist = [l.strip() for l in handle]
        training_genes[fold] = tmp_glist
    all_training = set([g
                        for fld in training_genes
                        for g in training_genes[fld]])
    return training_genes, all_training


def main(opts):
    # read iarc mutations
    df = pd.read_table(opts['input'])
    df['mutation'] = df['HGVSp_Short'].str[2:]

    # read preferred snvbox transcripts
    tx_df = pd.read_table(opts['preferred_tx'])

    # merge the two
    merged_df = pd.merge(df, tx_df, on='Hugo_Symbol', how='left')

    # add a UID column
    merged_df['UID'] = range(len(merged_df))

    # save output
    mycols = ['UID', 'transcript', 'mutation']
    if opts['cross_val_dir']:
        training_genes, all_training = read_cross_val(opts['cross_val_dir'])
        if not os.path.exists(opts['output']): os.mkdir(opts['output'])
        tmp = merged_df[~merged_df['Hugo_Symbol'].isin(all_training)]
        save_path = os.path.join(opts['output'], opts['prefix']+'_chasm_input_passenger.txt')
        tmp[mycols].to_csv(save_path, index=False, sep=' ', header=None)
        for i in range(10):
            test_genes = all_training - set(training_genes[i])
            tmp = merged_df[merged_df['Hugo_Symbol'].isin(test_genes)]
            save_path = os.path.join(opts['output'], opts['prefix']+'_chasm_input{0}.txt'.format(i))
            tmp[mycols].to_csv(save_path, index=False, header=None, sep=' ')
    else:
        merged_df[mycols].to_csv(opts['output'], sep=' ', index=False, header=None)


if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)
