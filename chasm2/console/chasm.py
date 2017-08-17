#!/usr/bin/env python
# fix problems with pythons terrible import system
import sys
import os
file_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(file_dir, '../'))
sys.path.append(os.path.join(file_dir, '../../'))

# package specific imports
import chasm2
import chasm2.python.utils as utils

# general imports
import argparse
import pandas as pd
import numpy as np
import logging

logger = logging.getLogger(__name__)  # module logger

def parse_arguments():
    # make a parser
    info = 'Predicts driver missense mutations in cancer'
    parent_parser = argparse.ArgumentParser(description=info)

    # logging arguments
    parent_parser.add_argument('-ll', '--log-level',
                               type=str,
                               action='store',
                               default='',
                               help='Write a log file (--log-level=DEBUG for debug mode, '
                               '--log-level=INFO for info mode)')
    parent_parser.add_argument('-l', '--log',
                               type=str,
                               action='store',
                               default='stdout',
                               help='Path to log file. (accepts "stdout")')
    parent_parser.add_argument('-v', '--verbose',
                               action='store_true',
                               default=False,
                               help='Flag for more verbose log output')

    # add subparsers
    subparsers = parent_parser.add_subparsers(title='Analysis types', dest='kind')
    help_info = 'Prepare snvbox input for snvGetGenomic'
    description = help_info + '. Separates training data into drivers vs passengers by a heuristic.'
    parser_prepSnvbox = subparsers.add_parser('prepSnvboxInput',
                                              help=help_info,
                                              description=description)
    help_info = 'Merge additional features into snvbox output'
    description = help_info + '.'
    parser_mergeFeat = subparsers.add_parser('mergeFeatures',
                                             help=help_info,
                                             description=description)
    help_info = 'Create null distribution for random forest scores'
    description = help_info + '.'
    parser_nullDist = subparsers.add_parser('nullDistribution',
                                            help=help_info,
                                            description=description)
    help_info = 'Calculates combined score using CHASM2 and 20/20+'
    description = help_info + '.'
    parser_combScore = subparsers.add_parser('combinedScore',
                                             help=help_info,
                                             description=description)
    help_info = 'Merge additional features into snvbox output for benchmark assessments'
    description = help_info + '.'
    parser_mergeBenchmarkFeat = subparsers.add_parser('mergeBenchmarkFeatures',
                                             help=help_info,
                                             description=description)

    # program arguments
    myparsers = [parser_prepSnvbox, parser_mergeFeat, parser_nullDist,
                 parser_combScore, parser_mergeBenchmarkFeat]
    for i, parser in enumerate(myparsers):
        # group of parameters
        major_parser = parser.add_argument_group(title='Major options')
        advance_parser = parser.add_argument_group(title='Advanced options')

        # set the CLI params
        if i==0:
            help_str = 'Somatic mutations in MAF format'
            major_parser.add_argument('-i', '--input',
                                    type=str, required=True,
                                    help=help_str)
            help_str = 'Directory containing results from mutsigcv v1.4'
            major_parser.add_argument('-m', '--mutsigcv-dir',
                                    type=str,
                                    help=help_str)
            help_str = 'File relating gene names to IDs'
            major_parser.add_argument('-g', '--gene-file',
                                    type=str, required=True,
                                    help=help_str)
            help_str = 'Snvbox file for driver mutations'
            major_parser.add_argument('-od', '--output-driver',
                                    type=str, required=True,
                                    help=help_str)
            help_str = 'Snvbox file for passenger mutations'
            major_parser.add_argument('-op', '--output-passenger',
                                    type=str, required=True,
                                    help=help_str)
        elif i==1:
            parser.add_argument('-i', '--input',
                                type=str, required=True,
                                help='Somatic mutations in MAF format')
            parser.add_argument('-ig', '--id2gene',
                                type=str, required=True,
                                help='file containing mapping of variant IDs to gene names')
            parser.add_argument('-s', '--snvbox',
                                type=str, required=True,
                                help='SNVbox features in tab delimited format')
            parser.add_argument('-hm', '--hotmaps',
                                type=str, required=True,
                                help='Hotmaps1d result files')
            parser.add_argument('-o', '--output',
                                type=str, required=True,
                                help='Output feature file with merged features')
        elif i==2:
            parser.add_argument('-t', '--ttplus-simdir',
                                type=str, required=True,
                                help='Directory containing simulation results for 20/20+')
            parser.add_argument('-c', '--chasm-simdir',
                                type=str, required=True,
                                help='Directory containing simulation results for CHASM2')
            parser.add_argument('-o', '--output',
                                type=str, required=True,
                                help='Null distribution file')
        elif i==3:
            parser.add_argument('-c', '--chasm2',
                                type=str, required=True,
                                help='Score results from CHASM2')
            parser.add_argument('-t', '--twentyTwentyPlus',
                                type=str, required=True,
                                help='Score results from 20/20+')
            parser.add_argument('-nd', '--null-distribution',
                                type=str, required=True,
                                help='Null distribution file')
            parser.add_argument('-o', '--output',
                                type=str, required=True,
                                help='Final output from CHASM2')
        elif i==4:
            parser.add_argument('-m', '--maf',
                                type=str, required=True,
                                help='Somatic mutations in MAF format')
            parser.add_argument('-s', '--snvbox',
                                type=str, required=True,
                                help='SNVbox features in tab delimited format')
            parser.add_argument('-b', '--benchmark',
                                type=str, required=True,
                                help='Benchmark mutations')
            parser.add_argument('-w', '--window',
                                type=str, required=True,
                                help='Window sizes for hotmaps')
            parser.add_argument('-n', '--null-distr-dir',
                                type=str, required=True,
                                help='Hotmaps1d null distribution')
            parser.add_argument('-o', '--output',
                                type=str, required=True,
                                help='Output feature file with merged features')
    args = parent_parser.parse_args()

    # handle logging
    if args.log_level or args.log:
        if args.log:
            log_file = args.log
        else:
            log_file = ''  # auto-name the log file
    else:
        log_file = os.devnull
    log_level = args.log_level
    utils.start_logging(log_file=log_file,
                        log_level=log_level,
                        verbose=args.verbose)  # start logging

    # create a dictionary for CLI options
    opts = vars(args)

    # log user entered command
    logger.info('Version: {0}'.format(chasm2.__version__))
    logger.info('Command: {0}'.format(' '.join(sys.argv)))
    return opts


def main(opts):
    if opts['kind'] == 'prepSnvboxInput':
        import chasm2.python.snvbox.make_snvbox_input as ms
        # add gene lists
        opts['oncogene'] = os.path.join(utils.data_dir, utils.config_dict['oncogene'])
        opts['tsg'] = os.path.join(utils.data_dir, utils.config_dict['tsg'])
        # run the snvbox creation code
        ms.main(opts)
    elif opts['kind'] == 'mergeFeatures':
        import chasm2.python.snvbox.add_gene_name as ag
        import chasm2.python.snvbox.merge_gene_features as gf
        import chasm2.python.snvbox.merge_hotmaps1d as hm
        # add gene name
        opts['features'] = None
        new_df = ag.main(opts)
        # add gene features
        opts['dgd'] = os.path.join(utils.data_dir, utils.config_dict['dgd'])
        opts['lawrence'] = os.path.join(utils.data_dir, utils.config_dict['mutsigcv'])
        opts['gene_length'] = os.path.join(utils.data_dir, utils.config_dict['gene_length'])
        opts['features'] = new_df
        new_df = gf.main(opts)
        # add hotmaps1d features
        opts['features'] = new_df
        hm.main(opts)
    elif opts['kind'] == 'mergeBenchmarkFeatures':
        import chasm2.python.snvbox.fetch_hotspot_pvalue as fp
        import chasm2.python.snvbox.merge_gene_features as gf
        import chasm2.python.snvbox.merge_hotmaps1d_benchmark as hm
        final_output = opts['output']
        snvbox = opts['snvbox']
        # fetch hotspot p-values
        opts['output'] = None
        opts['input'] = opts['benchmark']
        new_df = fp.main(opts)
        # add hotmaps1d features
        opts['hotmaps'] = new_df
        opts['features'] = opts['snvbox']
        new_df = hm.main(opts)
        # add gene features
        opts['dgd'] = os.path.join(utils.data_dir, utils.config_dict['dgd'])
        opts['lawrence'] = os.path.join(utils.data_dir, utils.config_dict['mutsigcv'])
        opts['gene_length'] = os.path.join(utils.data_dir, utils.config_dict['gene_length'])
        opts['features'] = new_df
        opts['output'] = final_output
        new_df = gf.main(opts)
    elif opts['kind'] == 'nullDistribution':
        import chasm2.python.null_dist as nd
        nd.main(opts)
    elif opts['kind'] == 'combinedScore':
        import chasm2.python.combined_score as cs
        cs.main(opts)


def cli_main():
    # run main with CLI options
    opts = parse_arguments()
    main(opts)


if __name__ == "__main__":
    cli_main()
