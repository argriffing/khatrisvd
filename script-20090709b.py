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
from khatrisvd import util

logging.basicConfig(level=logging.DEBUG)

def main():
    X = util.file_to_whitespace_separated_matrix('khatrisvd/fivetimes.txt')
    n = len(X)

    permutation = range(n)
    random.shuffle(permutation)
    X = heatmap.get_permuted_rows(X, permutation)

    root = treebuilder.build_tree(X)
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
