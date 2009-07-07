"""
Define some correlation related functions.

In this module, X is supposed to be a data matrix
with more rows than columns.
The motivation is that each row represents a
gene and each column represents a genetic line.
"""

import numpy as np

def khatri_rao_square(M):
    """
    @param M: a numpy array
    """
    np.vstack(np.kron(row, row) for row in M.T).T

def test():
    print 'ohai'
