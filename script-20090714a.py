"""
This is a junk script.

Test the box-filtering.
"""

import logging

import numpy as np

from khatrisvd import treebuilder
from khatrisvd import splitbuilder
from khatrisvd import khorr
from khatrisvd import heatmap
from khatrisvd import util
from khatrisvd import gradient

logging.basicConfig(level=logging.DEBUG)

def main():
    pathname_in = 'mmc-data-files/Starvation_Residual.TXT'
    X = util.file_to_comma_separated_matrix(pathname_in, has_headers=True)
    print X.shape
    # create the tree from the data
    root = treebuilder.build_tree(X)
    ordered_indices = root.ordered_labels()
    X = heatmap.get_permuted_rows(X, ordered_indices)
    # create the standardized data for drawing the small heatmap
    Z = khorr.get_standardized_matrix(X)
    color_function = gradient.correlation_to_rgb
    # show the big heatmap
    pathname_out = 'big.png'
    im = heatmap.get_heatmap_image(np.dot(Z, Z.T), color_function)
    fout = open(pathname_out, 'wb')
    im.save(fout)
    fout.close()
    # show the small heatmap
    pathname_out = 'small.png'
    im = heatmap.get_reduced_heatmap_image(Z, color_function, reduction=5)
    fout = open(pathname_out, 'wb')
    im.save(fout)
    fout.close()

if __name__ == '__main__':
    main()
