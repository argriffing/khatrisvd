"""
This is a junk script.

Maybe some stuff in this script can be salvaged for later use.
Read a data matrix from stdin.
Write a squared correlation cluster tree newick string to stdout.
The input matrix should have row and column headers,
and elements should be comma separated.
"""

import logging
import sys

import numpy as np

from khatrisvd import treebuilder

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
    logging.debug('creating the tree')
    root = treebuilder.build_tree(X)
    logging.debug('creating the newick string')
    print root.get_newick_string()

if __name__ == '__main__':
    main()

