"""
File: merge_benchmark.py
Author: Collin Tokheim
Email: ctokheim@jhu.edu
Github: ctokheim
Description: Merge the results of each of the methods
"""
import pandas as pd
import argparse
import os

def parse_arguments():
    info = 'Merge the results of each of the methods'
    parser = argparse.ArgumentParser(description=info)
    parser.add_argument('-b', '--benchmark-dir',
                        type=str, required=True,
                        help='Directory containing the benchmark data')
    args = parser.parse_args()
    return vars(args)


def read_candra(input_path, output_path):
    """Read candra results"""
    mycols = ['Chrom', 'Coordinate', 'Ref_Allele', 'Mut_Allele', 'Strand']
    candra_input_df = pd.read_table(input_path, header=None, names=mycols)
    candra_input_df['UID'] = range(len(candra_input_df))
    candra_input_df = candra_input_df.drop_duplicates(subset=mycols)
    candra_output_df = pd.read_table(output_path).drop_duplicates()
    merged = pd.merge(candra_input_df, candra_output_df, on=mycols, how='left')
    merged = merged.rename(columns={'CanDrA_Score': 'CanDrA'})
    out_cols = ['UID', 'CanDrA']
    return merged[out_cols]


def read_fathmm(input_path, output_path,):
    """Read FATHMM results."""
    fathmm_input_df = pd.read_table(input_path, header=None, sep=' ', names=['Protein ID', 'Substitution'])
    fathmm_output_df = pd.read_table(output_path)
    merged_fathmm = pd.merge(fathmm_input_df, fathmm_output_df.drop_duplicates(subset=['Protein ID', 'Substitution']),
                             on=['Protein ID', 'Substitution'], how='left')
    merged_fathmm = merged_fathmm.rename(columns={'Score': 'FATHMM'})
    return merged_fathmm


def read_chasm2(path):
    df = pd.read_table(path)
    df = df.rename(columns={'driver': 'CHASM2'})
    if 'CHASM2_genome' in df.columns:
        out_cols = ['UID', 'CHASM2', 'CHASM2_genome', 'CHASM2_pval',
                    'CHASM2_genome_pval', 'CHASM2_qval', 'CHASM2_genome_qval']
    else:
        out_cols = ['UID', 'CHASM2']
    return df[out_cols]


def read_chasm(path):
    """Read the CHASM results"""
    df = pd.read_table(path)
    df = df.rename(columns={'MutationID': 'UID', 'PValue':'CHASM_pval'})
    out_cols = ['UID', 'CHASM', 'CHASM_pval']
    return df[out_cols]


def read_annovar(path):
    """Read the annovar results."""
    df = pd.read_table(path, na_values=['.'])
    # rename columns
    rename_dict = {
        'VEST3_score': 'VEST',
        'CADD_raw': 'CADD',
        'Polyphen2_HDIV_score': 'Polyphen2_hdiv',
        'Polyphen2_HVAR_score': 'Polyphen2_hvar',
        'MutationAssessor_score': 'MutationAssessor',
        'SIFT_score': 'SIFT'
    }
    df = df.rename(columns=rename_dict)
    # add UID column
    df['UID'] = range(len(df))
    # output columns
    out_cols = ['UID', 'VEST', 'CADD', 'Polyphen2_hdiv', 'Polyphen2_hvar',
                'SIFT', 'MutationAssessor', 'SIFT', 'REVEL', 'MCAP']
    return df[out_cols]


def read_parssnp(path):
    df = pd.read_table(path, na_values=['.'])
    df['UID'] = range(len(df))
    out_cols = ['UID', 'ParsSNP']
    return df[out_cols]


def read_transfic(path):
    names = ['UID', 'gene', 'TransFIC', 'impact1', 'score2', 'impact2', 'score3', 'impact3']
    df = pd.read_table(path, header=None, names=names)
    out_cols = ['UID', 'TransFIC']
    return df[out_cols]


def read_input(benchmark_dir, benchmark):
    maf_like = ['msk_impact', 'patrick_et_al', 'mc3']
    if benchmark in maf_like:
        path = os.path.join(benchmark_dir, '{0}.maf'.format(benchmark))
        df = pd.read_table(path)
        df['UID'] = range(len(df))
        df = df.drop_duplicates(subset=['Hugo_Symbol', 'HGVSp_Short'])
        if 'Variant_Classification' in df.columns:
            df = df[df['Variant_Classification']=='Missense_Mutation']
    else:
        path = os.path.join(benchmark_dir, '{0}.txt'.format(benchmark))
        df = pd.read_table(path)
        df['UID'] = range(len(df))
    if benchmark == 'iarc_tp53':
        p53_re_cols = ['WAF1nWT', 'MDM2nWT', 'BAXnWT', 'h1433snWT', 'AIP1nWT',
                       'GADD45nWT', 'NOXAnWT', 'P53R2nWT']
        df['class'] = df[p53_re_cols].median(axis=1)

    return df


def main(opts):
    bench_dir = opts['benchmark_dir']
    maf_like = ['patrick_et_al', 'msk_impact', 'mc3']
    non_maf = ['berger_et_al', 'berger_et_al_egfr', 'kim_et_al', 'iarc_tp53']
    benchmarks = maf_like + non_maf

    for b in benchmarks:
        # read input file
        input_df = read_input(bench_dir, b)

        # read predictions
        candra_in = os.path.join(bench_dir, 'methods/input/{0}.candra_input.txt'.format(b))
        candra_out = os.path.join(bench_dir, 'methods/output/{0}.candra_output.txt'.format(b))
        candra_df = read_candra(candra_in, candra_out)
        candra_in = os.path.join(bench_dir, 'methods/input/{0}.candra_input.txt'.format(b))
        candra_out = os.path.join(bench_dir, 'methods/output/{0}.candra_plus_output.txt'.format(b))
        candra_plus_df = read_candra(candra_in, candra_out).rename(columns={'CanDrA': 'CanDrA plus'})
        chasm_path = os.path.join(bench_dir, 'methods/output/{0}.chasm_output.txt'.format(b))
        chasm_df = read_chasm(chasm_path)
        chasm2_path = os.path.join(bench_dir, 'methods/output/{0}.chasm2_output.txt'.format(b))
        if b != 'mc3': chasm2_df = read_chasm2(chasm2_path)
        transfic_path = os.path.join(bench_dir, 'methods/output/{0}.transfic_output.txt'.format(b))
        transfic_df = read_transfic(transfic_path)
        annovar_path = os.path.join(bench_dir, 'methods/output/{0}.annovar_output.hg19_multianno.txt'.format(b))
        annovar_df = read_annovar(annovar_path)
        parssnp_path = os.path.join(bench_dir, 'methods/output/ParsSNP.output.{0}.annovar_output.hg19_multianno.txt'.format(b))
        parssnp_df = read_parssnp(parssnp_path)
        annovar_df['ParsSNP'] = parssnp_df['ParsSNP']

        # merge input with annovar
        mycols = ['UID', 'Hugo_Symbol', "HGVSp_Short", 'class']
        merged_df = pd.merge(input_df[mycols], annovar_df,
                            on='UID', how='left')

        # merge results
        merged_df = pd.merge(merged_df, candra_df, on='UID', how='left')
        merged_df = pd.merge(merged_df, candra_plus_df, on='UID', how='left')
        merged_df = pd.merge(merged_df, chasm_df, on='UID', how='left')
        if b != 'mc3': merged_df = pd.merge(merged_df, chasm2_df, on='UID', how='left')
        merged_df = pd.merge(merged_df, transfic_df, on='UID', how='left')

        # save results
        output_path = os.path.join(bench_dir, 'methods/output/{0}_results.txt'.format(b))
        merged_df.to_csv(output_path, sep='\t', index=False)


if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)


