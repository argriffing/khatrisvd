"""
Draw a heatmap.

This uses PIL.
During testing it dumps an image in the current directory;
this might be a problem.
"""

import dendro

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


class HeatmapHelper:
    """
    This is a vehicle for the dendrogram line drawing function.
    """

    def __init__(self, image, dx, dy):
        """
        @param image: a PIL image
        @param dx: the dendrogram x offset
        @param dy: the dendrogram y offset
        """
        self.image = image
        self.dx = dx
        self.dy = dy
        #TODO add an option for the side where the dendrogram goes
        self.left = True

    def on_draw_dendrogram_line(self, line):
        """
        @param line: each endpoint is a (breadth_offset, height_offset) pair
        """
        #TODO draw a line using PIL instead of putting the pixels one by one
        ((a, b), (c, d)) = line
        if a == c:
            for height_offset in range(min(b, d), max(b, d) + 1):
                x = self.dx + height_offset
                y = self.dy + a
                self.image.putpixel((x, y), (0, 0, 0))
        elif b == d:
            for breadth_offset in range(min(a, c), max(a, c) + 1):
                x = self.dx + b
                y = self.dy + breadth_offset
                self.image.putpixel((x, y), (0, 0, 0))


def get_heatmap_with_dendrogram(M, root, filename):
    """
    The input matrix is meant to be an entrywise squared correlation matrix.
    The indices of M should have been ordered conformantly with the dendrogram order.
    @param M: a matrix with entries between zero and one
    @param root: the root of the tree used to create the dendrogram
    @param filename: where to put the image
    """
    #TODO add an option for the side where the dendrogram goes
    left = True
    #TODO draw in a way that doesn't suck
    n = len(M)
    # initialize the color blue
    blue_hue, blue_saturation, value = colorsys.rgb_to_hsv(0, 0, 1.0)
    saturation = blue_saturation
    # define dendrogram parameters
    breadth_gap = 2
    height_gap = 3
    dendrogram_breadth = dendro.get_dendrogram_breadth(root, breadth_gap)
    dendrogram_height = dendro.get_dendrogram_height(root, height_gap)
    dendrogram_dx = n*3 + 1
    dendrogram_dy = 1
    # define the width of the image
    image_width = n*3 + dendrogram_height + 2
    image_height = n*3
    # create the image using M as the saturation value
    myimage = Image.new('RGB', (image_width, image_height), 'white')
    # draw the dendrogram
    helper = HeatmapHelper(myimage, dendrogram_dx, dendrogram_dy)
    dendro.draw_dendrogram(root, breadth_gap, height_gap, helper.on_draw_dendrogram_line)
    # draw the heatmap
    for i in range(n):
        for j in range(n):
            saturation = M[i][j]
            colorsys_rgb = colorsys.hsv_to_rgb(blue_hue, saturation, value)
            pil_rgb = tuple(int(x*255) for x in colorsys_rgb)
            for u in range(3):
                for v in range(3):
                    x = i*3 + u
                    y = j*3 + v
                    myimage.putpixel((x, y), pil_rgb)
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
