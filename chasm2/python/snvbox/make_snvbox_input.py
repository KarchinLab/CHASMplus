"""
File: make_snvbox_input.py
Author: Collin Tokheim
Email: ctokheim@jhu.edu
Github: ctokheim
Description: Converts the MC3 MAF into input for snvbox
"""
import pandas as pd
import argparse
import os
import numpy as np
import csv


def parse_arguments():
    info = 'Converts the MC3 MAF into input for snvbox'
    parser = argparse.ArgumentParser(description=info)
    parser.add_argument('-i', '--input',
                        type=str, required=True,
                        help='maf file')
    parser.add_argument('-m', '--mutsigcv-dir',
                        type=str, default=None,
                        help='Directory containing driver predictions')
    parser.add_argument('-onco', '--oncogene',
                        type=str, required=True,
                        help='List containing oncogenes')
    parser.add_argument('-tsg', '--tsg',
                        type=str, required=True,
                        help='List containing tumor suppressors')
    parser.add_argument('-g', '--gene-file',
                        type=str, required=True,
                        help='File relating variant IDs to their gene')
    parser.add_argument('-od', '--output-driver',
                        type=str, required=True,
                        help='genomic SNVBox input format for drivers')
    parser.add_argument('-op', '--output-passenger',
                        type=str,
                        help='genomic SNVBox input format for passengers')
    args = parser.parse_args()
    return vars(args)


def read_signif_genes(base_dir, driver_list):
    result_dict = {}
    for myfile in os.listdir(base_dir):
        mypath = os.path.join(base_dir, myfile)
        cancer_type = myfile[:-4]
        if cancer_type == 'README': continue

        # get the driver list from result file
        tmp_df = pd.read_table(mypath, index_col=0)
        tmp_df_subset = tmp_df.loc[driver_list].copy()
        mydrivers = set(tmp_df_subset[tmp_df_subset['qvalue']<=.1].index.tolist())
        print('{0}: {1}'.format(cancer_type, len(mydrivers)))
        result_dict[cancer_type] = mydrivers
    return result_dict


def main(opts):
    # read file
    df = pd.read_table(opts['input'])

    # figure out high mutation samples
    if 'Tumor_Sample_Barcode' in df.columns:
        my_cts = df['Tumor_Sample_Barcode'].value_counts()
        high_count_samps = my_cts[my_cts>500].index.tolist()

    # read oncogenes/tsgs
    if opts['oncogene'] or opts['tsg']:
        with open(opts['oncogene']) as handle:
            ogs = [l.strip() for l in handle]
        with open(opts['tsg']) as handle:
            tsgs = [l.strip() for l in handle]
        driver = ogs+tsgs
    else:
        raise ValueError('Expected list of oncogenes and tumor suppressor genes')

    # process data
    df['mystrand'] = '+'
    #df['Start_Position'] = df['Start_Position'] - 1
    df['Chromosome'] = 'chr' + df['Chromosome'].astype(str)
    if 'ID' not in df.columns: df['ID'] = list(range(len(df)))
    if 'CODE' not in df.columns: df['CODE'] = ''

    # filter data frame
    cols = ['ID', 'Chromosome', 'Start_Position', #'End_Position',
            'mystrand', 'Reference_Allele', 'Tumor_Seq_Allele2'] #, 'oncogenic']
    is_single = (df['Reference_Allele'].str.len()==1) & (df['Tumor_Seq_Allele2'].str.len()==1)
    is_not_indel = (df['Reference_Allele']!='-') & (df['Tumor_Seq_Allele2']!='-')
    df = df[is_single & is_not_indel]
    df = df[df['Variant_Classification']=='Missense_Mutation']
    #df = df.drop_duplicates(['Transcript_ID', 'HGVSp_Short'])
    #df = df[~is_not_indel]

    # save snvbox input
    # filter to keep only mutations that are in genes
    # deemed significant for their corresponding cancer type
    if opts['mutsigcv_dir']:
        signif_gene_dict = read_signif_genes(opts['mutsigcv_dir'], driver)
        is_driver = []
        for i, row in df.iterrows():
            if row['Tumor_Sample_Barcode'] in high_count_samps:
                is_driver.append(False)
            elif row['Hugo_Symbol'] in signif_gene_dict.get(row['CODE'], []):
                is_driver.append(True)
            else:
                is_driver.append(False)
    else:
        is_driver = [True for i in range(len(df))]
    is_driver = pd.Series(is_driver, index=df.index)
    # just assume each gene is reported on one transcript if no transcript
    # column seen
    if 'Transcript_ID' not in df.columns: df['Transcript_ID'] = df['Hugo_Symbol']
    # save likely driver mutations
    df[is_driver].drop_duplicates(['Transcript_ID', 'HGVSp_Short'])[cols].to_csv(opts['output_driver'], sep='\t', index=False, header=False)
    # save likely passenger mutations
    driver_mut_list = set(map(tuple, df[is_driver].drop_duplicates(['Transcript_ID', 'HGVSp_Short'])[['Transcript_ID', 'HGVSp_Short']].values.tolist()))
    is_prev_driver = df.apply(lambda x, s: (x['Transcript_ID'], x['HGVSp_Short']) in s, args=(driver_mut_list,), axis=1)
    non_driver_df = df[(~is_driver) & (~is_prev_driver)].drop_duplicates(['Transcript_ID', 'HGVSp_Short'])
    non_driver_df[cols].to_csv(opts['output_passenger'], sep='\t', index=False, header=False)

    # save gene file
    cols2 = ['ID', 'Hugo_Symbol']
    df[cols2].to_csv(opts['gene_file'], sep='\t', index=False, header=False)


if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)
