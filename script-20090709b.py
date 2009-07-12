"""
This is a junk script.

Maybe it will be useful as an example.
"""

import random
import logging

import numpy as np

from khatrisvd import treebuilder
from khatrisvd import splitbuilder
from khatrisvd import khorr
from khatrisvd import heatmap

logging.basicConfig(level=logging.DEBUG)

def main():
    X = splitbuilder.get_data('khatrisvd/fivetimes.txt')
    n = len(X)

    permutation = range(n)
    random.shuffle(permutation)
    X = heatmap.get_permuted_rows(X, permutation)

    U, S = khorr.data_to_reduced_laplacian_sqrt(X)
    tree_data = treebuilder.TreeData()
    root = treebuilder.build_tree(U, S, range(len(U)), tree_data)
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
