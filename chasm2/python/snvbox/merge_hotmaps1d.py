"""
File: merge_hotmaps1d_feature.py
Author: Collin Tokheim
Email: ctokheim@jhu.edu
Github: ctokheim
Description: Merge the hotmaps1d features into snvbox features
"""
import pandas as pd
import argparse

# logging
import logging
logger = logging.getLogger(__name__)  # module logger


def parse_arguments():
    info = 'Merge the hotmaps1d features into snvbox features'
    parser = argparse.ArgumentParser(description=info)
    parser.add_argument('-i', '--input',
                        type=str, required=True,
                        help='Mutation file (MAF)')
    parser.add_argument('-hm', '--hotmaps',
                        type=str, required=True,
                        help='Hotmaps1d result file')
    parser.add_argument('-f', '--features',
                        type=argparse.FileType('r'), required=True,
                        help='Features file from snvget')
    parser.add_argument('-o', '--output',
                        type=str, required=True,
                        help='Feature file with merged hotmaps result')
    args = parser.parse_args()
    return vars(args)


def main(opts):
    # read in data
    #hotmaps_list = opts['hotmaps'].split(',')
    #width_list = [0, 5, 10]

    # handle mutations
    # figure out the column header
    with open(opts['input']) as handle:
        header = next(handle).split('\t')
        useful_cols = ['Gene', 'Protein_Change', 'Variant_Classification',
                       'Transcript_ID', 'Hugo_Symbol', 'Protein_position',
                       'HGVSp_Short', 'Tumor_Sample_Barcode']
        usecols = list(set(useful_cols) & set(header))
    mut_df = pd.read_table(opts['input'], usecols=usecols)
    # if not a full MAF fill in necessary info
    if 'Transcript_ID' not in mut_df.columns:
        mut_df = mut_df.rename(columns={'Gene': 'Hugo_Symbol',
                                        'Protein_Change': 'HGVSp_Short'})
        mut_df['Transcript_ID'] = mut_df['Hugo_Symbol']
        mut_df['Tumor_Sample_Barcode'] = '.'
        mut_df['Protein_position'] = mut_df['HGVSp_Short'].str.extract('^p.[A-Z]([0-9]+)[A-Z]$')
    # filter only missense
    mut_df['UID'] = range(len(mut_df))
    #mut_df2 = mut_df[mut_df['Variant_Classification']=='Missense_Mutation']

    # merge the two data frames
    logger.info('merging hotmaps. . .')
    hotmaps_df = pd.read_table(opts['hotmaps'])
    width_list = hotmaps_df['window length'].unique()
    for i, width in enumerate(width_list):
        #hotmaps_df = pd.read_table(hotmaps_list[i])
        tmp_hotmaps_df = hotmaps_df[hotmaps_df['window length']==width].copy()
        tmp_cols = ['Transcript_ID', 'Protein_position', 'Tumor_Sample_Barcode', 'UID']
        hmap_cols = ['index', 'p-value ({0})'.format(width)]
        tmp_hotmaps_df = tmp_hotmaps_df.rename(columns={'p-value': 'p-value ({0})'.format(width)})

        if i == 0:
            merged_df = pd.merge(tmp_hotmaps_df[hmap_cols], mut_df[tmp_cols],
                                 left_on='index', right_index=True, how='right')
        else:
            merged_df = pd.merge(tmp_hotmaps_df[hmap_cols], merged_df,
                                 on='index', how='right')

    # fill with mean of group
    logger.info('na filling . . .')
    cols_of_interest = ['Transcript_ID', 'Protein_position']
    grp = merged_df.groupby(cols_of_interest)
    pcols = ['p-value (0)', 'p-value (5)', 'p-value (10)']
    means = grp[pcols].mean().reset_index()
    merged_df = pd.merge(merged_df, means,
                         on=cols_of_interest, how='left',
                         suffixes=('', '_merged'))

    # fill in missing values
    for w in width_list:
        pcol_name = 'p-value ({0})'.format(w)
        # fill in with same codon
        isnull = merged_df[pcol_name].isnull()
        merged_df.loc[isnull, pcol_name] = merged_df.loc[isnull, pcol_name+'_merged']

        # fill in with overall mean if not found
        merged_df[pcol_name] = merged_df[pcol_name].fillna(merged_df[pcol_name].mean())
    #for w in width_list:
        #pcol_name = 'p-value ({0})'.format(w)
        #tmp_pval = grp[pcol_name].transform(lambda x: x.fillna(x.mean()))
        #merged_df[pcol_name] = tmp_pval

        # fill with mean overall if not found
        #merged_df[pcol_name] = merged_df[pcol_name].fillna(merged_df[pcol_name].mean())
    logger.info('done fill')

    # save class series
    logger.info('merge previous features ...')
    if opts['features'] is str:
        df = pd.read_table(opts['features'])
    else:
        df = opts['features']
    class_ser = df['class']
    del df['class']
    mycols = ['UID', #'gene',
              'p-value (0)', 'p-value (5)', 'p-value (10)']
    #final_df = pd.merge(df, merged_df.rename(columns={'index':'UID',
                                                      #'Hugo_Symbol': 'gene'
                                                      #})[mycols],
                        #on='UID', how='left')
    final_df = pd.merge(df, merged_df[mycols], on='UID', how='left')
    final_df['class'] = class_ser
    if opts['output']:
        final_df.to_csv(opts['output'], sep='\t', index=False)

    return final_df


if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)
