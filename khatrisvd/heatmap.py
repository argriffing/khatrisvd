"""
Draw a heatmap.

This uses PIL.
During testing it dumps an image in the current directory;
this might be a problem.
"""

import Image
import numpy as np

import colorsys
import unittest
import random

def get_inverse_permutation(order):
    """
    @param order: a permutation of range(n)
    @return: the inverse permutation
    """
    n = len(order)
    result = [0]*n
    for a, b in enumerate(order):
        result[b] = a
    return result

def get_permuted_rows(M, ordered_indices):
    """
    @param M: a matrix
    @param ordered_indices: maps final indices to original indices
    @return: a matrix whose rows have been reordered
    """
    return np.vstack([M[i] for i in ordered_indices])

def get_permuted_columns(M, ordered_indices):
    """
    @param M: a matrix
    @param ordered_indices: maps final indices to original indices
    @return: a matrix whose rows have been reordered
    """
    return get_permuted_rows(M.T, ordered_indices).T

def get_permuted_rows_and_columns(M, ordered_indices):
    """
    @param M: a matrix
    @param ordered_indices: maps final indices to original indices
    @return: a matrix whose rows and columns have been reordered
    """
    return get_permuted_rows(get_permuted_columns(M, ordered_indices), ordered_indices)

def get_heatmap(M, filename):
    """
    The input is meant to be an entrywise squared correlation matrix.
    The indices are expected to have been reordered already.
    @param M: a matrix with entries between zero and one
    @param filename: where to put the image
    """
    n = len(M)
    # initialize the color blue
    blue_hue, blue_saturation, value = colorsys.rgb_to_hsv(0, 0, 1.0)
    saturation = blue_saturation
    colorsys_rgb = colorsys.hsv_to_rgb(blue_hue, saturation, value)
    pil_rgb = tuple(int(x*255) for x in colorsys_rgb)
    # create the image using M as the saturation value
    myimage = Image.new('RGB', (n, n), 'red')
    for i in range(n):
        for j in range(n):
            saturation = M[i][j]
            colorsys_rgb = colorsys.hsv_to_rgb(blue_hue, saturation, value)
            pil_rgb = tuple(int(x*255) for x in colorsys_rgb)
            myimage.putpixel((i, j), pil_rgb)
    fout = open(filename, 'wb')
    myimage.save(fout)
    fout.close()


class TestMe(unittest.TestCase):

    def test_heatmap(self):
        """
        Make a test image.
        """
        filename = 'out.png'
        n = 100
        M = [[0]*n for i in range(n)]
        for i in range(20):
            for j in range(20):
                M[i][j] = .5
        for i in range(n):
            M[i][i] = 1
        get_heatmap(M, filename)

    def test_reordering(self):
        """
        Make two test images.
        """
        original_filename = 'out-original.png'
        reordered_filename = 'out-reordered.png'
        n = 100
        correlated_indices = set(random.sample(range(n), 20))
        # create the original heatmap
        M = [[0]*n for i in range(n)]
        for i in correlated_indices:
            for j in correlated_indices:
                M[i][j] = .5
        for i in range(n):
            M[i][i] = 1
        M = np.array(M)
        get_heatmap(M, original_filename)
        # create the reordered heatmap
        ordered_indices = list(correlated_indices) + list(set(range(n)) - correlated_indices)
        M_reordered = get_permuted_rows_and_columns(M, ordered_indices)
        get_heatmap(M_reordered, reordered_filename)

    def test_inverse_permutation(self):
        observed = get_inverse_permutation([2, 0, 1, 3])
        expected = [1, 2, 0, 3]
        self.assertEqual(expected, observed)


if __name__ == '__main__':
    unittest.main()
