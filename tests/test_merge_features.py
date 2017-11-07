# fix problems with pythons terrible import system
import os
import sys
file_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(file_dir, '..'))

import chasm2.console.chasm as chasm

def test_merge_features():
    """Merge in additional features"""
    # prepare dictionary of parameters
    maf_path = os.path.join(file_dir, 'data/mutations.hg38.maf')
    id2gene_path = os.path.join(file_dir, 'result/id2gene.txt')
    hotmaps_path = os.path.join(file_dir, 'result/hotmaps1d/result.txt')
    snvbox_feature_path = os.path.join(file_dir, 'result/snvbox_features.txt')
    output_path = os.path.join(file_dir, 'result/snvbox_features_merged.txt')
    opts = {
        'kind': 'mergeFeatures',
        'input': maf_path,
        'id2gene': id2gene_path,
        'hotmaps': hotmaps_path,
        'snvbox': snvbox_feature_path,
        'output': output_path,
    }

    # run prep snvbox code
    chasm.main(opts)
