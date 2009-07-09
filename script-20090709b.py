"""
This is a junk script.

Maybe it will be useful as an example.
"""

import numpy as np
import random

from khatrisvd import treebuilder
from khatrisvd import splitbuilder
from khatrisvd import khorr
from khatrisvd import heatmap

def main():
    X = splitbuilder.get_data('khatrisvd/fivetimes.txt')
    n = len(X)

    permutation = range(n)
    random.shuffle(permutation)
    X = heatmap.get_permuted_rows(X, permutation)

    L_sqrt = khorr.data_to_laplacian_sqrt(X)
    tree_data = treebuilder.TreeData(splitbuilder.split_svd, treebuilder.update_svd)
    root = treebuilder.build_tree(L_sqrt, range(len(X)), tree_data)
    # show the unordered heatmap
    filename = 'unordered.png'
    RoR = np.corrcoef(X)**2
    heatmap.get_heatmap(RoR, filename)
    # show the ordered heatmap
    filename = 'reordered.png'
    ordered_indices = root.ordered_labels()
    M = heatmap.get_permuted_rows_and_columns(RoR, ordered_indices)
    heatmap.get_heatmap(M, filename)

if __name__ == '__main__':
    main()
