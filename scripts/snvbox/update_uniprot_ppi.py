"""
File: update_exon_feat.py
Author: Collin Tokheim
Email: ctokheim@jhu.edu
Github: ctokheim
Description: Add ppi feature to uniprot_features table
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
        uniprot_ix, pos_ix, count_ix = header.index('uniprot_id'), header.index('residue'), header.index('count')

        # run mysql command to reset to zero
        myquery = (
        'UPDATE Uniprot_features '
        '    SET insider_ppi=0 '
        )
        cursor.execute(myquery)
        # commit changes
        db.commit()

        # iterate over each line of the file
        for line in myreader:
            # get the relevant columns
            uniprot, pos, count = line[uniprot_ix], line[pos_ix], line[count_ix]

            # run mysql command
            myquery = (
            'UPDATE Uniprot_features '
            '    SET insider_ppi={count} '
            'WHERE Acc=\'{uniprot_id}\' AND Pos={pos}'
            ).format(count=count, uniprot_id=uniprot, pos=int(pos))
            cursor.execute(myquery)
            # commit changes
            db.commit()

            num_rows_affected = cursor.rowcount
            if num_rows_affected == 0:
                colnames = (
                    'Acc',  'Pos', 'BINDING', 'ACT_SITE', 'SITE',
                    'LIPID', 'METAL', 'CARBOHYD', 'DNA_BIND', 'NP_BIND',
                    'CA_BIND', 'DISULFID', 'SE_CYS', 'MOD_RES', 'PROPEP',
                    'SIGNALP', 'TRANSMEM', 'COMPBIAS', 'REP',  'MOTIF',
                    'ZN_FING', 'REGIONS', 'PPI',  'RNABD', 'TF',
                    'LOC',  'MMBRBD', 'Chrom', 'PostModRec',
                    'PostModEnz', 'insider_ppi'
                )
                # run mysql command
                myquery = (
                'INSERT INTO Uniprot_features ({colnames}) '
                'VALUES (\'{uniprot_id}\', {pos}, {zeros}, {count})'
                ).format(colnames=', '.join(colnames), uniprot_id=uniprot,
                         pos=int(pos), zeros=', '.join(['0']*(len(colnames)-3)),
                         count=count)
                cursor.execute(myquery)
                # commit changes
                db.commit()



if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)
