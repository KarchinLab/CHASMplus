"""
File: analyze_rarity.py
Author: Collin Tokheim
Email: ctokheim@jhu.edu
Github: ctokheim
Description: Analyze the rarity of the significant mutations
"""
import pandas as pd
import argparse


def parse_arguments():
    info = 'Analyze the rarity of the significant mutations'
    parser = argparse.ArgumentParser(description=info)
    parser.add_argument('-c', '--chang-et-al',
                        type=str, required=True,
                        help='Chang et al result')
    parser.add_argument('-m', '--maf',
                        type=str, required=True,
                        help='Mutation file')
    parser.add_argument('-o', '--output',
                        type=str, required=True,
                        help='Output file')
    args = parser.parse_args()
    return vars(args)


def main(opts):
    # read data
    rename_dict = {'Amino_Acid_Position': 'Protein_position'}
    chang_df = pd.read_table(opts['chang_et_al']).rename(columns=rename_dict)
    mut_df = pd.read_table(opts['maf'])

    # count the number of samples
    num_samps = mut_df.groupby('CODE')['Tumor_Sample_Barcode'].nunique()

    # figure out which mutations are significant
    mut_cols = ['Hugo_Symbol', 'Transcript_ID', 'Protein_position', 'HGVSp_Short']
    chang_df = pd.merge(chang_df, mut_df[mut_cols],
                        on=['Hugo_Symbol', 'Protein_position'], how='left')
    signif = chang_df[chang_df['qvalue']<.01]
    signif_list = (signif['Hugo_Symbol'] + '_' + signif['Protein_position'].astype(str)).tolist()

    # filter the mutations
    mut_label = (mut_df['Hugo_Symbol'] + '_' + mut_df['Protein_position'].astype(str))
    mut_df['mutation'] = mut_label
    prot_change = (mut_df['Hugo_Symbol'] + '_' + mut_df['Protein_position'].astype(str))
    mut_df['prot_change'] = prot_change
    signif_mut_df = mut_df[mut_df['prot_change'].isin(signif_list)]

    # calc the fraction of samples
    count = signif_mut_df.groupby(['mutation', 'CODE'])['Tumor_Sample_Barcode'].nunique()
    frac = count / num_samps
    frac_df = frac.reset_index(name='fraction of samples')

    # calculate the maximum fraction
    pivot = frac_df.pivot(index='mutation', columns='CODE', values='fraction of samples')
    max_frac = pivot.max(axis=1)
    max_cancer_type = pivot.idxmax(axis=1)
    tmp_dict = {'fraction of samples': max_frac, 'cancer type': max_cancer_type}
    max_df = pd.DataFrame(tmp_dict).reset_index()

    # merge in the total mutation count
    mut_ct = signif_mut_df['mutation'].value_counts()
    mut_ct = mut_ct.reset_index(name='number of mutations').rename(columns={'index': 'mutation'})
    max_df = pd.merge(max_df, mut_ct, on='mutation', how='left')

    # label the different categories
    is_rare = (max_df['fraction of samples']<0.01) | (max_df['number of mutations']==1)
    is_common = max_df['fraction of samples']>.05
    max_df['category'] = 'intermediate'
    max_df.loc[is_common, 'category'] = 'common'
    max_df.loc[is_rare, 'category'] = 'rare'

    # save output
    max_df.to_csv(opts['output'], sep='\t', index=False)


if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)


