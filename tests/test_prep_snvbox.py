# fix problems with pythons terrible import system
import os
import sys
file_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(file_dir, '..'))

import chasm2.console.chasm as chasm

def test_prep_snvbox():
    """Test the preparation of SNVBox input."""
    # prepare dictionary of parameters
    maf_path = os.path.join(file_dir, 'data/mutations.hg38.maf')
    id2gene_path = os.path.join(file_dir, 'result/id2gene.txt')
    driver_path = os.path.join(file_dir, 'result/driver.snvbox_input.txt')
    passenger_path = os.path.join(file_dir, 'result/passenger.snvbox_input.txt')
    opts = {
        'kind': 'prepSnvboxInput',
        'input': maf_path,
        'gene_file': id2gene_path,
        'output_driver': driver_path,
        'output_passenger': passenger_path,
        'mutsigcv_dir': None
    }

    # run prep snvbox code
    chasm.main(opts)
