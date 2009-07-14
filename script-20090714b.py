"""
This is a junk script.

Run the program on the big matrix.
python script.py <data> <newick>

"""

import logging
import sys

import numpy as np

from khatrisvd import treebuilder
from khatrisvd import splitbuilder
from khatrisvd import khorr
from khatrisvd import heatmap
from khatrisvd import util
from khatrisvd import gradient
from khatrisvd import progress

logging.basicConfig(level=logging.DEBUG)

def main():
    # first read the data file
    pathname_in = sys.argv[1]
    X = util.file_to_comma_separated_matrix(pathname_in, has_headers=True)
    logging.debug('matrix size: ' + str(X.shape))
    # now read the newick file
    pathname_in = sys.argv[2]
    fin = open(pathname_in)
    newick_string = fin.read()
    fin.close()
    ordered_indices = heatmap.get_newick_order(newick_string)
    logging.debug('newick order size: ' + str(len(ordered_indices)))
    # reorder the rows of the data file
    X = heatmap.get_permuted_rows(X, ordered_indices)
    # create the standardized data for drawing the small heatmap
    Z = khorr.get_standardized_matrix(X)
    color_function = gradient.correlation_to_rgb
    # initialize the progress bar
    reduction = 100
    p = progress.Bar(len(X)/reduction+1)
    # show the small heatmap
    pathname_out = 'small.png'
    im = heatmap.get_reduced_heatmap_image(Z, color_function, reduction=reduction, line_callback=p.update)
    fout = open(pathname_out, 'wb')
    im.save(fout)
    fout.close()
    # finish the progress bar
    p.finish()

if __name__ == '__main__':
    main()
