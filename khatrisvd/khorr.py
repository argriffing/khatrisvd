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
import unittest
import math

import numpy as np

import util

def get_fiedler_vector(U, S):
    """
    The first element of the output vector is forced to be nonnegative.
    @param U: columns are eigenvectors
    @param S: conformant singular values
    @return: the eigenvector corresponding to the smallest nonsmall eigenvalue
    """
    best_w, best_i = min((w, i) for i, w in enumerate(S) if not util.is_small(w))
    v = U.T[best_i]
    if v[0] < 0:
        return -v
    else:
        return v

def remove_small_vectors(U, S):
    """
    Remove negligible values in S and corresponding columns in U.
    The inputs are not directly modified, but the results are returned.
    @param U: columns are eigenvectors
    @param S: conformant singular values
    @return: U_reduced, S_reduced
    """
    columns = [U.T[i] for i, w in enumerate(S) if not util.is_small(w)]
    U_reduced = np.vstack(columns).T
    S_reduced = [w for w in S if not util.is_small(w)]
    return U_reduced, S_reduced

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
    p, n = M.shape
    sqrt2 = math.sqrt(2)
    R_p = p
    R_n = (n*(n+1))/2
    R = np.zeros((R_p, R_n))
    for row_index, row in enumerate(M):
        index = 0
        R_row = R[row_index]
        for i, a in enumerate(row):
            R_row[index] = a*a
            index += 1
            for b in row[i+1:]:
                R_row[index] = sqrt2 * a * b
                index += 1
    return R

def get_standardized_matrix(X):
    """
    Get the square root of the correlation matrix, in a certain sense.
    If X is a data matrix and Z is the returned matrix,
    then ZZ' should be the correlation matrix of rows of X.
    @param M: a data matrix
    @return: a standardized matrix where each row has mean 0 and var 1/ncols
    """
    nrows, ncols = X.shape
    Z = util.get_row_centered_matrix(X)
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
    logging.debug('standardized_to_augmented_C: doing a singular value decomposition')
    U, S_array, VT = np.linalg.svd(Z, full_matrices=0)
    Z_reduced = np.array([row[:-1] for row in U*S_array])
    logging.debug('standardized_to_augmented_C: creating the khatri rao row square')
    W = reduced_khatri_rao_row_square(Z_reduced)
    return W

def data_to_laplacian_sqrt(X):
    """
    If the output is U, S then (U*S)(U*S)' is like a Laplacian.
    @param X: a data matrix
    @return: U, S
    """
    logging.debug('data_to_laplacian_sqrt: creating the standardized matrix')
    Z = get_standardized_matrix(X)
    logging.debug('data_to_laplacian_sqrt: creating the augmented matrix')
    Q = standardized_to_augmented_C(Z)
    logging.debug('data_to_laplacian_sqrt: creating the column centered matrix')
    W = util.get_column_centered_matrix(Q)
    logging.debug('data_to_laplacian_sqrt: manually cleaning up old matrices')
    del Z
    del Q
    logging.debug('data_to_laplacian_sqrt: doing a singular value decomposition')
    U, S_array, VT = np.linalg.svd(W, full_matrices=0)
    S_pinv_array = np.array([0 if util.is_small(x) else 1/x for x in S_array])
    return U, S_pinv_array

def data_to_reduced_laplacian_sqrt(X):
    """
    @param X: a data matrix
    @return: a matrix whose product with its transpose is like a Laplacian
    """
    U, S = data_to_laplacian_sqrt(X)
    return remove_small_vectors(U, S)


class TestMe(unittest.TestCase):

    def assertAllClose(self, A, B, msg=None):
        """
        Assert that two matrices are almost entrywise identical.
        @param A: a matrix
        @param B: a matrix
        @param msg: shown when the assertion fails
        """
        self.assertTrue(np.allclose(A, B), msg=msg)

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
        expected_RoR = util.hadamard_square(np.corrcoef(X))
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
    
    def test_fiedler(self):
        p = 100
        n = 7
        X = np.random.random((p, n))
        # get a matrix that we are treating like the laplacian
        HDH = util.get_doubly_centered_matrix(np.corrcoef(X) ** 2)
        L = np.linalg.pinv(HDH)
        U, S, VT = np.linalg.svd(L, full_matrices=0)
        expected_fiedler = get_fiedler_vector(U, S)
        # get the square root of the summed laplacian without using pxp operations
        U, S = data_to_laplacian_sqrt(X)
        self.assertAllClose(np.dot(U*S, (U*S).T), L)
        observed_fiedler = get_fiedler_vector(U, S)
        self.assertAllClose(expected_fiedler, observed_fiedler)
        # get another square root of the summed laplacian
        U, S = remove_small_vectors(U, S)
        self.assertAllClose(np.dot(U*S, (U*S).T), L)
        self.assertAllClose(expected_fiedler, observed_fiedler)
        # get yet another square root of the summed laplacian
        U, S, VT = np.linalg.svd(U*S, full_matrices=0)
        self.assertAllClose(np.dot(U*S, (U*S).T), L)
        self.assertAllClose(expected_fiedler, observed_fiedler)

    def test_laplacian_summation_shortcut(self):
        p = 100
        n = 7
        k = 3
        X = np.random.random((p, n))
        # get a matrix that we are treating like the laplacian
        HDH = util.get_doubly_centered_matrix(np.corrcoef(X) ** 2)
        L = np.linalg.pinv(HDH)
        L_summed = util.sum_last_rows_and_columns(L, k)
        U, S, VT = np.linalg.svd(L_summed, full_matrices=0)
        expected_fiedler = get_fiedler_vector(U, S)
        # get the square root of the summed laplacian without using pxp operations
        U, S = data_to_laplacian_sqrt(X)
        B = util.sum_last_rows(U*S, k)
        self.assertAllClose(np.dot(B, B.T), L_summed)
        # get another square root
        U, S = remove_small_vectors(U, S)
        B = util.sum_last_rows(U*S, k)
        U, S, VT = np.linalg.svd(B, full_matrices=0)
        observed_fiedler = get_fiedler_vector(U, S)
        self.assertAllClose(np.dot(U*S, (U*S).T), L_summed)
        self.assertAllClose(expected_fiedler, observed_fiedler)

    def test_extended_laplacian_shortcut(self):
        """
        Demo commutativity of SVD and summation under certain conditions.
        """
        p = 100
        n = 7
        k = 3
        # define sets of arbitrary row indices
        first_rows = set((3, 5, 9))
        second_rows = set((12, 13, 2))
        # make a random data matrix
        X = np.random.random((p, n))
        # get a matrix that we are treating like the laplacian
        HDH = util.get_doubly_centered_matrix(np.corrcoef(X) ** 2)
        L = np.linalg.pinv(HDH)
        L_summed = util.sum_arbitrary_rows_and_columns(L, first_rows)
        L_summed = util.sum_arbitrary_rows_and_columns(L_summed, second_rows)
        U, S, VT = np.linalg.svd(L_summed, full_matrices=0)
        expected_fiedler = get_fiedler_vector(U, S)
        for f in (data_to_laplacian_sqrt, data_to_reduced_laplacian_sqrt):
            # get the square root of the laplacian without using pxp operations
            U, S = f(X)
            A = U*S
            # sum rows of the summed laplacian sqrt
            B = util.sum_arbitrary_rows(A, first_rows)
            # get the first criterion matrix
            U, S_array, VT = np.linalg.svd(B, full_matrices=0)
            QB = U*S_array
            # sum rows of the first criterion matrix
            QB_summed = util.sum_arbitrary_rows(QB, second_rows)
            # sum more rows of the summed laplacian sqrt
            C = util.sum_arbitrary_rows(B, second_rows)
            # get the second criterion matrix directly from the summed laplacian sqrt
            U, S_array, VT = np.linalg.svd(C, full_matrices=0)
            QC_direct_fiedler = get_fiedler_vector(U, S_array)
            QC_direct = U*S_array
            # get the second criterion matrix from the summed first criterion matrix
            U, S_array, VT = np.linalg.svd(QB_summed, full_matrices=0)
            QC_indirect_fiedler = get_fiedler_vector(U, S_array)
            QC_indirect = U*S_array
            # check the equivalence of the matrix squares
            self.assertAllClose(np.dot(QC_direct, QC_direct.T), L_summed)
            self.assertAllClose(np.dot(QC_indirect, QC_indirect.T), L_summed)
            self.assertAllClose(expected_fiedler, QC_direct_fiedler)
            self.assertAllClose(expected_fiedler, QC_indirect_fiedler)

if __name__ == '__main__':
    unittest.main()
