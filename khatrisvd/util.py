"""
Various simple matrix functions.
"""

import unittest

import numpy as np


def is_small(x, eps=1e-12):
    return abs(x) < eps

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

def hadamard_square(M):
    """
    @param M: a matrix
    @return: the Hadamard square of M
    """
    return M*M

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

def sum_last_rows_and_columns(M, k):
    """
    Get the matrix that is M with its last k rows and columns summed.
    @param M: a matrix
    @param k: the number of trailing rows and columns to sum
    @return: the modified matrix
    """
    return sum_last_rows(sum_last_columns(M, k), k)

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

def sum_arbitrary_rows_and_columns(M, indices_to_sum):
    """
    Get the matrix that is M with some rows and columns summed and moved to the end.
    The order of the columns that are not summed will be stable in the resulting matrix.
    @param M: a matrix
    @param column_indices_to_sum: the set of indices of rows and columns to sum
    """
    return sum_arbitrary_rows(sum_arbitrary_columns(M, indices_to_sum), indices_to_sum)


class TestMe(unittest.TestCase):

    def assertAllClose(self, A, B, msg=None):
        """
        Assert that two matrices are almost entrywise identical.
        @param A: a matrix
        @param B: a matrix
        @param msg: shown when the assertion fails
        """
        self.assertTrue(np.allclose(A, B), msg=msg)

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
