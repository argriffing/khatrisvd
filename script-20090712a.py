"""
This is for testing changes to the reduced khatri rao squaring function.
"""

import logging
import sys
import time

import numpy as np

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
    start_time = time.time()
    filename = sys.argv[1]
    logging.debug('reading input lines')
    fin = open(filename)
    lines = fin.readlines()
    fin.close()
    logging.debug('converting input lines to data matrix')
    X = read_matrix(lines)
    logging.debug('standardizing the data matrix')
    Z = khorr.get_standardized_matrix(X)
    logging.debug('augmenting the data matrix')
    W = khorr.standardized_to_augmented_C(Z)
    end_time = time.time()
    print end_time - start_time
    raw_input('kbye\n')

if __name__ == '__main__':
    main()

