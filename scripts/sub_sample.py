"""
File: sub_sample.py
Author: Collin Tokheim
Email: ctokheim@jhu.edu
Github: ctokheim
Description: Sub-sample the mutations
"""
import pandas as pd
import numpy as np
import os
import argparse


def parse_arguments():
    info = 'Sub-sample the mutations'
    parser = argparse.ArgumentParser(description=info)
    parser.add_argument('-i', '--input',
                        type=str, required=True,
                        help='MAF file')
    parser.add_argument('-n', '--num-samps',
                        type=int, required=True,
                        help='Number of times to sub-sample')
    parser.add_argument('-s', '--seed',
                        type=int, default=101,
                        help='Seed of pseudo random number generator (Default: 101)')
    parser.add_argument('-o', '--output',
                        type=str, required=True,
                        help='Output directory')
    args = parser.parse_args()
    return vars(args)


def main(opts):
    # read mutations
    df = pd.read_table(opts['input'])
    num_samps = opts['num_samps']

    # sub-sample
    unique_samps = df.Tumor_Sample_Barcode.unique()
    num_unique = len(unique_samps)
    prng = np.random.RandomState(seed=opts['seed'])
    probs = prng.uniform(size=num_samps)
    for i, prob in enumerate(probs):
        # figure out the samples
        samp_size = int(np.ceil(prob * num_unique))
        tmp_samples = prng.choice(unique_samps, size=samp_size, replace=False)

        # create directory
        dir_path = os.path.join(opts['output'], 'sample{0}'.format(i))
        if not os.path.exists(dir_path): os.makedirs(dir_path)

        # save mutations
        save_path = os.path.join(dir_path, 'mutations.maf')
        df[df.Tumor_Sample_Barcode.isin(tmp_samples)].to_csv(save_path, sep='\t', index=False)

    # save the sample rates
    samp_rate_df = pd.DataFrame({'sample': range(len(probs)),
                                 'fraction': probs,
                                 'total samples': num_unique})
    save_path = os.path.join(opts['output'], 'sample_rate.txt')
    samp_rate_df.to_csv(save_path, sep='\t', index=False)



if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)

