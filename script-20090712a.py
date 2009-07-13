"""
This is for testing changes to the reduced khatri rao squaring function.
"""

import logging
import sys
import time

import numpy as np

from khatrisvd import khorr
from khatrisvd import util

logging.basicConfig(level=logging.DEBUG)

def main():
    start_time = time.time()
    filename = sys.argv[1]
    logging.debug('reading input lines')
    fin = open(filename)
    lines = fin.readlines()
    fin.close()
    logging.debug('converting input lines to data matrix')
    X = util.lines_to_comma_separated_matrix(lines, has_headers=True)
    logging.debug('standardizing the data matrix')
    Z = khorr.get_standardized_matrix(X)
    logging.debug('augmenting the data matrix')
    W = khorr.standardized_to_augmented_C(Z)
    end_time = time.time()
    print 'nseconds:', end_time - start_time
    print 'W.size:', W.size
    print 'W.itemsize:', W.itemsize
    print 'W.size * W.itemsize:', W.size * W.itemsize
    raw_input('kbye\n')

if __name__ == '__main__':
    main()

