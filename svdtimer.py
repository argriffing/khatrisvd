"""
This is just a script to test the svd speed.

Reasonable inputs are n=450 p=10000.
"""

import time
import sys

import numpy as np

def main():
    usage = 'Usage: ' + sys.argv[0] + ' <n> <p>'
    if len(sys.argv) != 3:
        print usage
        return
    script_name, nrows_string, ncols_string = sys.argv
    n, p = int(nrows_string), int(ncols_string)
    print 'create a', n, 'by', p, 'matrix:'
    start_time = time.time()
    X = np.random.random((n, p))
    print time.time() - start_time
    print 'take the svd of the matrix:'
    start_time = time.time()
    U, S, VT = np.linalg.svd(X, full_matrices=0)
    print time.time() - start_time

if __name__ == '__main__':
    main()
