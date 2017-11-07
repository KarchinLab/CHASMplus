# fix problems with pythons terrible import system
import os
import sys
file_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(file_dir, '..'))

import chasm2.console.chasm as chasm

def test_combined_score():
    """Test the combined score code."""
    # prepare dictionary of parameters
    chasm2_path = os.path.join(file_dir, 'result/chasm2_result_pretrained.txt')
    ttplus_path = os.path.join(file_dir, 'result/pretrained_output/results/r_random_forest_prediction.txt')
    null_path = os.path.join(file_dir, 'result/chasm2_null_distribution_pretrained.txt')
    output_path = os.path.join(file_dir, 'result/chasm2_result_pretrained_final.txt')
    opts = {
        'kind': 'combinedScore',
        'chasm2': chasm2_path,
        'twentyTwentyPlus': ttplus_path,
        'null_distribution': null_path,
        'output': output_path,
    }

    # run prep snvbox code
    chasm.main(opts)
