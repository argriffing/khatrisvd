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

import logging

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
    Z = get_row_centered_matrix(X)
    Z = np.vstack(row / np.std(row) for row in Z)
    Z /= math.sqrt(ncols)
    return Z

def get_row_centered_matrix(M):
    """
    @param M: a matrix
    @return: each row sums to zero
    """
    return np.vstack((row - np.mean(row)) for row in M)

def get_column_centered_matrix(M):
    """
    @param M: a matrix
    @return: each column sums to zero
    """
    return get_row_centered_matrix(M.T).T

def get_doubly_centered_matrix(M):
    """
    @param M: a matrix
    @return: each row and column sums to zero
    """
    return get_row_centered_matrix(get_column_centered_matrix(M))

def sum_last_rows(M, k):
    """
    Get the matrix that is M with its last k rows summed.
    @param M: a matrix
    @param k: the number of trailing rows to sum
    @return: the modified matrix
    """
    return np.vstack([M[:-k], np.sum(M[-k:], axis=0)])

def sum_last_columns(M, k):
    """
    Get the matrix that is M with its last k columns summed.
    @param M: a matrix
    @param k: the number of trailing columns to sum
    @return: the modified matrix
    """
    return sum_last_rows(M.T, k).T

def sum_arbitrary_rows(M, row_indices_to_sum):
    """
    Get the matrix that is M with some rows summed and moved to the end.
    The last row in the output matrix will be equal to the sum of the specified rows.
    The order of the rows that are not summed will be stable in the resulting matrix.
    @param M: a matrix
    @param row_indices_to_sum: the set of indices of rows to sum
    """
    n = len(M)
    row_indices_to_keep = set(range(n)) - set(row_indices_to_sum)
    row_sum = sum(M[i] for i in row_indices_to_sum)
    new_rows = [M[i] for i in sorted(row_indices_to_keep)] + [row_sum]
    return np.vstack(new_rows)

def sum_arbitrary_columns(M, column_indices_to_sum):
    """
    Get the matrix that is M with some columns summed and moved to the end.
    The last column in the output matrix will be equal to the sum of the specified columns.
    The order of the columns that are not summed will be stable in the resulting matrix.
    @param M: a matrix
    @param column_indices_to_sum: the set of indices of columns to sum
    """
    return sum_arbitrary_rows(M.T, column_indices_to_sum).T

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
    logging.debug('standardized_to_augmented_C: doing a singular value decomposition')
    U, S_array, VT = np.linalg.svd(Z, full_matrices=0)
    Z_reduced = np.array([row[:-1] for row in U*S_array])
    logging.debug('standardized_to_augmented_C: creating the khatri rao row square')
    W = reduced_khatri_rao_row_square(Z_reduced)
    return W

def is_small(x, eps=1e-12):
    return abs(x) < eps

def data_to_laplacian_sqrt(X):
    """
    @param X: a data matrix
    @return: a matrix whose product with its transpose is like a Laplacian
    """
    logging.debug('data_to_laplacian_sqrt: creating the standardized matrix')
    Z = get_standardized_matrix(X)
    logging.debug('data_to_laplacian_sqrt: creating the augmented matrix')
    Q = standardized_to_augmented_C(Z)
    logging.debug('data_to_laplacian_sqrt: creating the column centered matrix')
    W = get_column_centered_matrix(Q)
    logging.debug('data_to_laplacian_sqrt: manually cleaning up old matrices')
    del Z
    del Q
    logging.debug('data_to_laplacian_sqrt: doing a singular value decomposition')
    U, S_array, VT = np.linalg.svd(W, full_matrices=0)
    S_pinv_array = np.array([0 if is_small(x) else 1/x for x in S_array])
    L_sqrt = U*S_pinv_array
    return L_sqrt


class TestMe(unittest.TestCase):

    def assertAllClose(self, A, B, msg=None):
        """
        Assert that two matrices are almost entrywise identical.
        @param A: a matrix
        @param B: a matrix
        @param msg: shown when the assertion fails
        """
        self.assertTrue(np.allclose(A, B), msg=msg)

    def test_sum_arbitrary_rows(self):
        M = np.array([[1,2,3],[4,5,6],[7,8,9]])
        # first summation
        expected = np.array([[1,2,3],[11,13,15]])
        observed = sum_arbitrary_rows(M, (1, 2))
        self.assertAllClose(observed, expected)
        # second summation
        expected = np.array([[7,8,9],[5,7,9]])
        observed = sum_arbitrary_rows(M, (0, 1))
        self.assertAllClose(observed, expected)

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

    def test_trailing_row_and_column_summation(self):
        M = np.array([[1,2,3],[4,5,6],[7,8,9]])
        # test trailing row summation
        expected = np.array([[1,2,3],[11,13,15]])
        observed = sum_last_rows(M, 2)
        self.assertAllClose(observed, expected)
        # test trailing column summation
        expected = np.array([[1,5],[4,11],[7,17]])
        observed = sum_last_columns(M, 2)
        self.assertAllClose(observed, expected)

    def test_laplacian_summation_shortcut(self):
        p = 100
        n = 7
        k = 3
        X = np.random.random((p, n))
        # get a matrix that we are treating like the laplacian
        HDH = get_doubly_centered_matrix(np.corrcoef(X) ** 2)
        L = np.linalg.pinv(HDH)
        L_summed = sum_last_rows(sum_last_columns(L, k), k)
        # get the square root of the summed laplacian without using pxp operations
        L_sqrt = data_to_laplacian_sqrt(X)
        L_summed_sqrt = sum_last_rows(L_sqrt, k)
        # assert that the square of the square root is the summed laplacian
        L_summed_reconstructed = np.dot(L_summed_sqrt, L_summed_sqrt.T)
        self.assertAllClose(L_summed_reconstructed, L_summed)

    def test_extended_laplacian_shortcut(self):
        """
        Demo commutativity of SVD and summation under certain conditions.
        """
        p = 100
        n = 7
        k = 3
        # define sets of arbitrary rows
        first_rows = set((3, 5, 9))
        second_rows = set((12, 13, 2))
        # make a random data matrix
        X = np.random.random((p, n))
        # get a matrix that we are treating like the laplacian
        HDH = get_doubly_centered_matrix(np.corrcoef(X) ** 2)
        L = np.linalg.pinv(HDH)
        L_summed = sum_arbitrary_rows(sum_arbitrary_columns(L, first_rows), first_rows)
        L_summed = sum_arbitrary_rows(sum_arbitrary_columns(L_summed, second_rows), second_rows)
        # get the square root of the laplacian without using pxp operations
        A = data_to_laplacian_sqrt(X)
        # sum rows of the summed laplacian sqrt
        B = sum_arbitrary_rows(A, first_rows)
        # get the first criterion matrix
        U, S_array, VT = np.linalg.svd(B, full_matrices=0)
        QB = U*S_array
        # sum rows of the first criterion matrix
        QB_summed = sum_arbitrary_rows(QB, second_rows)
        # sum more rows of the summed laplacian sqrt
        C = sum_arbitrary_rows(B, second_rows)
        # get the second criterion matrix directly from the summed laplacian sqrt
        U, S_array, VT = np.linalg.svd(C, full_matrices=0)
        QC_direct = U*S_array
        # get the second criterion matrix from the summed first criterion matrix
        U, S_array, VT = np.linalg.svd(QB_summed, full_matrices=0)
        QC_indirect = U*S_array
        # check the equivalence
        self.assertAllClose(np.dot(QC_direct, QC_direct.T), L_summed)
        self.assertAllClose(np.dot(QC_indirect, QC_indirect.T), L_summed)

if __name__ == '__main__':
    unittest.main()
