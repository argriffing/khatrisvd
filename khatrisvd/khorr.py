"""
Define some functions related to correlation matrices.

In this module, X is supposed to be a data matrix
with more rows than columns.
The motivation is that each row represents a
gene and each column represents a genetic line.
Whenever a matrix is mentioned here,
it is represented as a two dimensional numpy array
unless otherwise specified.
"""

import numpy as np

import unittest
import math

def khatri_rao_square(M):
    """
    Compute the Khatri-Rao product of a matrix with itself
    @param M: a matrix
    @return: the Khatri-Rao product of M with itself
    """
    return khatri_rao_row_square(M.T).T

def khatri_rao_row_square(M):
    """
    This function is to rows as the Khatri-Rao square is to columns.
    @param M: a matrix
    @return: a matrix larger than M
    """
    return np.vstack(np.kron(r, r) for r in M)

def reduced_khatri_rao_row_square(M):
    """
    This is an adjustment of the Khatri-Rao row square.
    The redundant columns in the output are combined.
    @param M: a matrix
    @return: a matrix larger than M
    """
    augmented_rows = []
    for row in M:
        augmented_row = []
        for i, a in enumerate(row):
            augmented_row.append(a*a)
            for b in row[i+1:]:
                augmented_row.append(math.sqrt(2) * a * b)
        augmented_rows.append(np.array(augmented_row))
    return np.vstack(augmented_rows)

def hadamard_square(M):
    """
    @param M: a matrix
    @return: the Hadamard square of M
    """
    return M*M

def get_standardized_matrix(X):
    """
    Get the square root of the correlation matrix, in a certain sense.
    If X is a data matrix and Z is the returned matrix,
    then ZZ' should be the correlation matrix of rows of X.
    @param M: a data matrix
    @return: a standardized matrix where each row has mean 0 and var 1/ncols
    """
    nrows, ncols = X.shape
    Z = np.vstack((row - np.mean(row)) for row in X)
    Z = np.vstack(row / np.std(row) for row in Z)
    Z /= math.sqrt(ncols)
    return Z

def standardized_to_augmented_A(Z):
    """
    This is the first of three similar functions.
    Each of these functions returns a square root,
    in a certain sense, of the Hadamard squared correlation matrix.
    If Z has p rows and n columns,
    then the returned matrix will have p rows and n*n columns.
    @param Z: a standardized data matrix
    @return: an augmented matrix
    """
    return khatri_rao_row_square(Z)

def standardized_to_augmented_B(Z):
    """
    This is the second of three similar functions.
    Each of these functions returns a square root,
    in a certain sense, of the Hadamard squared correlation matrix.
    If Z has p rows and n columns,
    then the returned matrix will have p rows and n*(n+1)/2 columns.
    This function corrects for the redundant columns
    returned by the first function.
    @param Z: a standardized data matrix
    @return: an augmented matrix
    """
    return reduced_khatri_rao_row_square(Z)

def standardized_to_augmented_C(Z):
    """
    This is the third of three similar functions.
    Each of these functions returns a square root,
    in a certain sense, of the Hadamard squared correlation matrix.
    If Z has p rows and n columns,
    then the returned matrix will have p rows and n*(n-1)/2 columns.
    This function corrects for the rank reduction induced by standardization.
    @param Z: a standardized data matrix
    @return: an augmented matrix
    """
    U, S_array, VT = np.linalg.svd(Z, full_matrices=0)
    Z_reduced = np.array([row[:-1] for row in U*S_array])
    return reduced_khatri_rao_row_square(Z_reduced)


class TestMe(unittest.TestCase):

    def assertAllClose(self, A, B, msg=None):
        """
        Assert that two matrices are almost entrywise identical.
        @param A: a matrix
        @param B: a matrix
        @param msg: shown when the assertion fails
        """
        self.assertTrue(np.allclose(A, B), msg=msg)

    def test_hadamard_square(self):
        M = np.array([[1,2],[3,4]])
        expected = np.array([[1,4],[9,16]])
        observed = hadamard_square(M)
        self.assertAllClose(observed, expected)

    def test_khatri_rao_square(self):
        M = np.array([[1,2,3],[4,5,6]])
        expected = np.array([
            [ 1,  4,  9],
            [ 4, 10, 18],
            [ 4, 10, 18],
            [16, 25, 36]])
        observed = khatri_rao_square(M)
        self.assertAllClose(observed, expected)

    def test_correlation(self):
        X = np.random.random((20, 3))
        Z = get_standardized_matrix(X)
        observed = np.dot(Z, Z.T)
        expected = np.corrcoef(X)
        self.assertAllClose(observed, expected)

    def test_standardization_rank_reduction(self):
        """
        Demonstrate that standardization reduces the rank of a matrix.
        Also show how this affects the singular value decomposition.
        """
        X = np.random.random((20, 3))
        Z = get_standardized_matrix(X)
        U, S_array, VT = np.linalg.svd(Z, full_matrices=0)
        a, b, c = S_array
        self.assertTrue(a > b)
        self.assertTrue(b > c)
        self.assertAlmostEqual(c, 0)

    def test_standardized_to_augmented(self):
        """
        Show 3 ways to get a sqrt of a Hadamard squared correlation matrix.
        """
        p = 20
        n = 3
        X = np.random.random((p, n))
        Z = get_standardized_matrix(X)
        expected_RoR = hadamard_square(np.corrcoef(X))
        functions_and_ncols = [
                (standardized_to_augmented_A, n*n),
                (standardized_to_augmented_B, (n*(n+1))/2),
                (standardized_to_augmented_C, (n*(n-1))/2)]
        for f, expected_ncols in functions_and_ncols:
            W = f(Z)
            observed_nrows, observed_ncols = W.shape
            self.assertEqual(observed_ncols, expected_ncols, msg=f.__name__)
            observed_RoR = np.dot(W, W.T)
            self.assertAllClose(observed_RoR, expected_RoR, msg=f.__name__)


if __name__ == '__main__':
    unittest.main()
