"""
File: snvbox2fathmm.py
Author: Collin Tokheim
Email: ctokheim@jhu.edu
Github: ctokheim
Description: Generate protein input format for fathmm
"""
import pandas as pd
import MySQLdb
import argparse
import csv


def parse_arguments():
    info = 'Generate protein input format for fathmm'
    parser = argparse.ArgumentParser(description=info)
    parser.add_argument('-i', '--input',
                        type=str, required=True,
                        help='SNVbox transcript input format')
    parser.add_argument('--from-features',
                        action='store_true', default=False,
                        help='Input is feature file from snvget output')
    parser.add_argument('-mh', '--mysql-host',
                        type=str,
                        default='karchin-db01.icm.jhu.edu',
                        help='Host name for mysql')
    parser.add_argument('-mdb', '--mysql-db',
                        type=str,
                        default='SNVBox_20161028_sandbox',
                        help='Database name for snvbox')
    parser.add_argument('--mysql-user',
                        type=str, required=True,
                        help='MySQL user name')
    parser.add_argument('--mysql-passwd',
                        type=str, required=True,
                        help='MySQL password')
    parser.add_argument('-o', '--output',
                        type=argparse.FileType('w'), required=True,
                        help='FATHMM input format')
    args = parser.parse_args()
    return vars(args)


def main(opts):
    # open db connection
    db = MySQLdb.connect(host=opts['mysql_host'],
                         user=opts['mysql_user'],
                         passwd=opts['mysql_passwd'],
                         db=opts['mysql_db'])
    cursor = db.cursor()

    # read in input file
    cols = ['UID', 'transcript', 'mutation']
    if not opts['from_features']:
        df = pd.read_table(opts['input'], sep=' ', header=None, names=cols)
    else:
        df = pd.read_table(opts['input'], sep='\t')
        mut = df['ID'].str.extract('_([A-Z0-9]+)$')
        tx = df['ID'].str.extract('([A-Z0-9_.]+)_[A-Z0-9]+$')
        df['transcript'] = tx
        df['mutation'] = mut
        df = df[cols].copy()

    # iterate through each row
    output_list = []
    for ix, row in df.iterrows():
        # figure out which column to use in mysql table
        if row['transcript'].startswith('ENST'):
            col_name = 'EnsT'
        elif row['transcript'].startswith('NM_'):
            col_name = 'RefseqT'

        # construct mysql query
        myquery = "SELECT * FROM Transcript WHERE {0}='{1}'".format(col_name, row['transcript'])
        cursor.execute(myquery)

        # get uniprot ID
        tmp_result = cursor.fetchone()
        uniprot_id = tmp_result[-2]

        output_list.append([uniprot_id, row['mutation']])

    # save output
    with opts['output'] as handle:
        mywriter = csv.writer(handle, delimiter=' ', lineterminator='\n')
        mywriter.writerows(output_list)
    db.close()


if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)

