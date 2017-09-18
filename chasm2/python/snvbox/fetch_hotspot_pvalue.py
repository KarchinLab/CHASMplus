import pandas as pd
import argparse
import os


def parse_arguments():
    info = 'Retrieve the hotspot p-value from the null distribution from a set of queried mutations'
    parser = argparse.ArgumentParser(description=info)
    parser.add_argument('-i', '--input',
                        type=str, required=True,
                        help='Mutations to get hotspot p-value for')
    parser.add_argument('-m', '--maf',
                        type=str, required=True,
                        help='mutations in maf format')
    parser.add_argument('-w', '--window',
                        type=str, required=True,
                        help='Window size to consider')
    parser.add_argument('-n', '--null-distr-dir',
                        type=str, required=True,
                        help='Directory containing null distribution files')
    parser.add_argument('-o', '--output',
                        type=str, required=True,
                        help='result file with p-value')
    args = parser.parse_args()
    return vars(args)


def count2pval(mycount, count_list, pval_list):
    # return 0 p-value if higher than observed
    if mycount > count_list[0]: return 0

    # figure out the p-value
    for i in range(len(count_list)):
        if mycount == count_list[i]:
            return pval_list[i]
        if mycount > count_list[i]:
            return pval_list[i]

    # if 0 counts, return p-value of 1
    return 1.0


def main(opts):
    # set window size for hotspot
    window = list(map(int, opts['window'].split(',')))
    # set directory for null dist
    null_dir = opts['null_distr_dir']

    # read data
    df = pd.read_table(opts['input'])
    rename_dict = {
        'gene': 'Hugo_Symbol',
        'mutation' : 'HGVSp_Short',
        'ProtDescription': 'HGVSp_Short'
    }
    df = df.rename(columns=rename_dict)

    # add index
    if 'index' not in df.columns:
        df['index'] = range(len(df))

    # filtering
    if 'Variant_Classification' in df.columns:
        #num_aas = df.HGVSp_Short.str.findall('[A-Za-z]').str.len()
        is_missense_format = df.HGVSp_Short.str.match('^p.[A-Z][0-9]+[A-Z]$')
        is_missense = df['Variant_Classification']=='Missense_Mutation'
        #df = df[is_missense & (num_aas==2)]
        df = df[is_missense & is_missense_format]

    df['Protein_position'] = df['HGVSp_Short'].str[3:-1].astype(int)
    query_genes = df['Hugo_Symbol'].unique()

    # read maf
    maf_df = pd.read_table(opts['maf'])
    # keep only mutations in relevant genes
    maf_df = maf_df[maf_df['Hugo_Symbol'].isin(query_genes)].copy()
    # keep only missense
    is_valid_missense = (maf_df['Variant_Classification']=='Missense_Mutation') & \
                        (maf_df['Protein_position'].str.contains('^[0-9]+$'))
    maf_df = maf_df[is_valid_missense].copy()
    maf_df['Protein_position'] = maf_df['Protein_position'].astype(int)

    # go over each row
    mut_dict = {w: [] for w in window}
    pval_dict = {w: [] for w in window}
    missing_genes = []
    for ix, row in df.iterrows():
        # figure out the mutation counts
        is_gene = maf_df['Hugo_Symbol']==row['Hugo_Symbol']
        maf_df_gene = maf_df[is_gene]

        for w in window:
            is_window = (maf_df['Protein_position'] <= (row['Protein_position']+w)) & \
                        (maf_df['Protein_position'] >= (row['Protein_position']-w))
            num_mut = len(maf_df_gene[is_window])
            mut_dict[w].append(num_mut)

            # figure out the p-value
            mygene = row['Hugo_Symbol']
            null_path = os.path.join(null_dir, mygene+'.{0}.txt'.format(w))
            if os.path.exists(null_path):
                null_df = pd.read_table(null_path)
                tmp_ct_list, tmp_pval_list = null_df['mutation_count'].tolist(), null_df['p-value'].tolist()
                pval = count2pval(num_mut, tmp_ct_list, tmp_pval_list)
                pval_dict[w].append(pval)
            else:
                pval_dict[w].append(1)
                if mygene not in missing_genes:
                    missing_genes.append(mygene)

    # add columns
    for w in window:
        df['window_count ({0})'.format(w)] = mut_dict[w]
        df['p-value ({0})'.format(w)] = pval_dict[w]

    # save file
    if opts['output']:
        df.to_csv(opts['output'], sep='\t', index=False)

    #print('Missing: {0}'.format(len(missing_genes)))
    return df


if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)
