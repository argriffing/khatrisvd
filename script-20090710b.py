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
from khatrisvd import util

logging.basicConfig(level=logging.DEBUG)

def main():
    # read the data matrix from stdin
    logging.debug('reading input lines')
    lines = sys.stdin.readlines()
    logging.debug('converting input lines to data matrix')
    X = util.lines_to_comma_separated_matrix(lines, has_headers=True)
    logging.debug('creating the tree')
    root = treebuilder.build_tree(X)
    logging.debug('creating the newick string')
    print root.get_newick_string()

if __name__ == '__main__':
    main()

