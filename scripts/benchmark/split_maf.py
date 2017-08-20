"""
File: split_maf.py
Author: Collin Tokheim
Email: ctokheim@jhu.edu
Github: ctokheim
Description: Split the MAF file based on the training folds of CHASM
"""
import pandas as pd
import argparse
import glob
import os


def parse_arguments():
    info = 'Split the MAF file based on the training folds of CHASM'
    parser = argparse.ArgumentParser(description=info)
    parser.add_argument('-m', '--maf',
                        type=str, required=True,
                        help='MAF file to split')
    parser.add_argument('-c', '--cross-val-dir',
                        type=str, required=True,
                        help='Directory containing the genes used in each fold')
    parser.add_argument('-o', '--output',
                        type=str, required=True,
                        help='Output directory to save splits')
    args = parser.parse_args()
    return vars(args)


def main(opts):
    # read mutation data
    mut_df = pd.read_table(opts['maf'])
    mut_df['ID'] = range(len(mut_df))

    # read training folds
    pattern = os.path.join(opts['cross_val_dir'], 'training*.txt')
    training_genes = dict()
    for f in glob.glob(pattern):
        fold = int(f[-5])
        with open(f) as handle:
            tmp_glist = [l.strip() for l in handle]
        training_genes[fold] = tmp_glist
    all_training = set([g
                        for fld in training_genes
                        for g in training_genes[fld]])

    # go over each fold
    if not os.path.exists(opts['output']): os.mkdir(opts['output'])
    tmp = mut_df[~mut_df['Hugo_Symbol'].isin(all_training)]
    save_path = opts['output'] + '_split_passenger.txt'
    tmp.to_csv(save_path, index=False, sep='\t')
    for i in range(10):
        test_genes = all_training - set(training_genes[i])
        tmp = mut_df[mut_df['Hugo_Symbol'].isin(test_genes)]
        save_path = opts['output'] + '_split_driver{0}.txt'.format(i)
        tmp.to_csv(save_path, sep='\t', index=False)


if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)

