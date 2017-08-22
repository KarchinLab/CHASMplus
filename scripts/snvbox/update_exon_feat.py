"""
File: update_exon_feat.py
Author: Collin Tokheim
Email: ctokheim@jhu.edu
Github: ctokheim
Description: Add uniprot feature to Exon_Features table
"""
import MySQLdb
import csv
import argparse


def parse_arguments():
    info = 'Add uniprot features to Exon_Features table'
    parser = argparse.ArgumentParser(description=info)
    parser.add_argument('-i', '--input',
                        type=str, required=True,
                        help='Text file containing uniprot density feature')
    parser.add_argument('-mh', '--mysql-host',
                        type=str,
                        default='karchin-db01.icm.jhu.edu',
                        help='Host name for mysql')
    parser.add_argument('-mdb', '--mysql-db',
                        type=str,
                        default='SNVBox_dev_sandbox',
                        help='Database name for snvbox')
    parser.add_argument('--mysql-user',
                        type=str, required=True,
                        help='MySQL user name')
    parser.add_argument('--mysql-passwd',
                        type=str, required=True,
                        help='MySQL password')
    args = parser.parse_args()
    return vars(args)


def main(opts):
    # connect to mysql
    db = MySQLdb.connect(host=opts['mysql_host'],
                         user=opts['mysql_user'],
                         passwd=opts['mysql_passwd'],
                         db=opts['mysql_db'])
    cursor = db.cursor()

    with open(opts['input']) as handle:
        myreader = csv.reader(handle, delimiter='\t')
        header = next(myreader)
        uid_ix, exon_ix, den_ix = header.index('UID'), header.index('Exon'), header.index('UniprotDensity')

        # iterate over each line of the file
        for line in myreader:
            # get the relevant columns
            uid, exon, den = line[uid_ix], line[exon_ix], line[den_ix]

            # run mysql command
            myquery = (
            'UPDATE Exon_Features '
            '    SET uniprot_den={den} '
            'WHERE UID={uid} AND Exon={exon}'
            ).format(den=den, uid=uid, exon=exon)
            cursor.execute(myquery)

        # commit changes
        db.commit()


if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)
