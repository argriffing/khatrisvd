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

def file_to_comma_separated_matrix(filename, has_headers=False):
    fin = open(filename)
    lines = fin.readlines()
    fin.close()
    return _parse_lines(lines, has_headers=has_headers, use_comma=True)

def file_to_whitespace_separated_matrix(filename, has_headers=False):
    fin = open(filename)
    lines = fin.readlines()
    fin.close()
    return _parse_lines(lines, has_headers=has_headers, use_comma=False)

def lines_to_comma_separated_matrix(lines, has_headers=False):
    return _parse_lines(lines, has_headers=has_headers, use_comma=True)

def lines_to_whitespace_separated_matrix(lines, has_headers=False):
    return _parse_lines(lines, has_headers=has_headers, use_comma=False)

def _parse_lines(lines, has_headers=False, use_comma=False):
    """
    Turn raw lines of input into a numpy data array.
    Each row of data is supposed to represent a gene or a DNA tile.
    Each column represents an individual or genetic line or experimental condition.
    Each entry is an expression level or log ratio to a reference expression level or something.
    Of course, this same procedure could be applied to different data types.
    Entries, including headers, should be separated by commas or whitespace or something.
    @param lines: raw lines of input
    @param use_comma: True if comma separated; False if whitespace separated
    @param has_headers: True iff rows and columns have headers
    @return: a numpy array
    """
    # read the input file
    lines = [line.strip() for line in lines]
    # skip empty lines
    lines = [line for line in lines if line]
    # possibly skip the first line
    if has_headers:
        lines = lines[1:]
    # get rows of elements
    if use_comma:
        rows = [line.split(',') for line in lines]
    else:
        rows = [line.split() for line in lines]
    # possibly skip the first element of each line
    if has_headers:
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
