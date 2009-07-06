"""
This is just a script to test the svd speed.
"""

import numpy
import time

def main():
    n = 450;
    p = 10000;
    print 'create a', n, 'by', p, 'matrix:'
    start_time = time.time()
    X = numpy.random.random((n, p))
    print time.time() - start_time
    print 'take the svd of the matrix:'
    start_time = time.time()
    U, S, VT = numpy.linalg.svd(X, full_matrices=0)
    print time.time() - start_time

if __name__ == '__main__':
    main()
