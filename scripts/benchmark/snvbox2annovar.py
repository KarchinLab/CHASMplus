"""
File: snvbox2annovar.py
Author: Collin Tokheim
Email: ctokheim@jhu.edu
Github: ctokheim
Description: Script to convert SNVBox genomic format to annovar
"""
import pandas as pd
import argparse


def parse_arguments():
    info = 'Script to convert SNVBox genomic format to annovar'
    parser = argparse.ArgumentParser(description=info)
    parser.add_argument('-i', '--input',
                        type=str, required=True,
                        help='SNVBox genomic format file')
    parser.add_argument('-o', '--output',
                        type=str, required=True,
                        help='Annovar format file')
    args = parser.parse_args()
    return vars(args)


def comp(letters):
    complement = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    tmp = []
    for let in letters:
        tmp.append(complement[let])
    return ''.join(tmp)


def main(opts):
    # read data
    cols = ['UID', 'chrom', 'pos', 'strand', 'ref', 'alt']
    df = pd.read_table(opts['input'], header=None, names=cols)

    # get positive strand if negative
    df.loc[df['strand']=='-', 'ref'] = df.loc[df['strand']=='-', 'ref'].apply(comp)
    df.loc[df['strand']=='-', 'alt'] = df.loc[df['strand']=='-', 'alt'].apply(comp)

    # get end position
    df['end'] = df['pos'] + df['ref'].str.len() - 1

    # fix chrom names
    df['chrom_new'] = df['chrom'].str[3:]

    # save file
    mycols = ['chrom_new', 'pos', 'end', 'ref', 'alt']
    df[mycols].to_csv(opts['output'], sep='\t', index=False, header=None)


if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)


