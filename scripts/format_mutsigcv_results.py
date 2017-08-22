"""
File: format_mutsigcv_results.py
Author: Collin Tokheim
Email: ctokheim@jhu.edu
Github: ctokheim
Description: Format the MutSigCV results into standard cancer driver format
"""
import os
import pandas as pd
import glob
import argparse


def parse_arguments():
    info = 'Fromat the mutsigcv results into standard cancer driver format'
    parser = argparse.ArgumentParser(description=info)
    parser.add_argument('-i', '--input',
                        type=str, required=True,
                        help='Directory containing raw mutsigcv results')
    parser.add_argument('-o', '--output',
                        type=str, required=True,
                        help='Directory containing formatted output')
    args = parser.parse_args()
    return vars(args)


def main(opts):
    rename_dict = {
        'p': 'pvalue', 'q': 'qvalue'
    }
    out_cols = ['gene', 'transcript', 'protein_change',
                'score', 'pvalue', 'qvalue', 'info']
    for f in glob.glob(opts['input']+'/*'):
        cancer_type = os.path.basename(f)
        tmp_path = os.path.join(f, 'mutations.sig_genes.txt')
        df = pd.read_table(tmp_path)
        df = df.rename(columns=rename_dict)
        df['transcript'] = '.'
        df['protein_change'] = '.'
        df['score'] = '.'
        df['info'] = '.'

        # save output
        if not os.path.exists(opts['output']): os.makedirs(opts['output'])
        out_path = os.path.join(opts['output'], cancer_type+'.txt')
        df[out_cols].to_csv(out_path, sep='\t', index=False)


if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)

