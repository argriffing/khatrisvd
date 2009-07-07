"""
Define some correlation related functions.

In this module, X is supposed to be a data matrix
with more rows than columns.
The motivation is that each row represents a
gene and each column represents a genetic line.
"""

import numpy as np

import unittest

def khatri_rao_square(M):
    """
    @param M: a numpy array
    @return: the Khatri-Rao product of the matrix with itself
    """
    return np.vstack(np.kron(row, row) for row in M.T).T


class TestMe(unittest.TestCase):

    def test_khatri_rao_square(self):
        M = np.array([[1,2,3],[4,5,6]])
        expected = np.array([
            [ 1,  4,  9],
            [ 4, 10, 18],
            [ 4, 10, 18],
            [16, 25, 36]])
        observed = khatri_rao_square(M)
        self.assertTrue(np.allclose(observed, expected))


if __name__ == '__main__':
    unittest.main()
