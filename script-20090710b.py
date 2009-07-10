"""
This is a junk script.

Maybe some stuff in this script can be salvaged for later use.
Read a squared correlation cluster tree from stdin
and write a newick tree to stdout.
"""

import logging
import optparse
import sys

import numpy as np

from khatrisvd import treebuilder
from khatrisvd import splitbuilder
from khatrisvd import khorr

logging.basicConfig(level=logging.DEBUG)

def read_matrix(lines):
    """
    The input file is csv with headers.
    @param lines: raw lines of input
    """
    # read the input file
    lines = [line.strip() for line in lines]
    # skip empty lines
    lines = [line for line in lines if line]
    # skip the first line
    lines = lines[1:]
    # get rows of elements
    rows = [line.split(',') for line in lines]
    # skip the first element of each line
    rows = [row[1:] for row in rows]
    # convert elements to floats
    X = []
    for row_index, row in enumerate(rows):
        try:
            float_row = [float(x) for x in row]
        except ValueError, e:
            message_lines = [
                    'invalid number on data row %d' % (row_index+1),
                    str(e)]
            raise ValueError('\n'.join(message_lines))
    rows = [[float(x) for x in row] for row in rows]
    # return the matrix
    return np.array(rows)

def main():
    # read the data matrix from stdin
    logging.debug('reading input lines')
    lines = sys.stdin.readlines()
    logging.debug('converting input lines to data matrix')
    X = read_matrix(lines)
    # make the first split
    logging.debug('creating the sqrt laplacian matrix')
    L_sqrt = khorr.data_to_laplacian_sqrt(X)
    # make subsequent splits
    logging.debug('creating the tree')
    tree_data = treebuilder.TreeData(splitbuilder.split_svd, treebuilder.update_svd)
    root = treebuilder.build_tree(L_sqrt, range(len(X)), tree_data)
    # write the newick tree to stdout
    logging.debug('creating the newick string')
    print root.get_newick_string()

if __name__ == '__main__':
    main()

