"""
Draw a heatmap.

This uses PIL.
During testing it dumps an image in the current directory;
this might be a problem.
"""

import colorsys
import unittest
import random

import Image
import numpy as np

import dendro
import mtree
import gradient

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

def get_reduced_heatmap_image(M, entry_to_rgb, reduction=5, line_callback=None):
    """
    @param M: a tall matrix such that dot products of rows are mapped to colors
    @param entry_to_rgb: a function that returns an RGB color given an entry of MM'
    @param reduction: a square this many pixels across is reduced to a single pixel
    @param line_callback: a function that is called with the number of lines completed
    """
    # int((((x+1)/2)*63)+1/2)
    # int(x*31.5 + 32)
    alpha = 31.5
    beta = 32
    # precalculate the 64 rgb colors from the gradient data
    colors = [np.array([int(x*255) for x in rgb]) for rgb in gradient.g_correlation_gradient]
    # get the number of rows in the input matrix
    n = len(M)
    # get the number of pixels on one side of the output image
    npixels, remainder = divmod(n, reduction)
    if remainder:
        npixels += 1
    # draw the image pixel by pixel
    im = Image.new('RGB', (npixels, npixels), 'red')
    for xpixel_index in range(npixels):
        xblock = M[xpixel_index*reduction : (xpixel_index+1)*reduction]
        for ypixel_index in range(npixels):
            if xpixel_index <= ypixel_index:
                yblock = M[ypixel_index*reduction : (ypixel_index+1)*reduction]
                submatrix = np.array(alpha*np.dot(xblock, yblock.T) + beta, dtype=np.int32)
                rgb_sum = sum(colors[i] for i in submatrix.flat)
                rgb = tuple(rgb_sum / submatrix.size)
                im.putpixel((xpixel_index, ypixel_index), rgb)
                im.putpixel((ypixel_index, xpixel_index), rgb)
        if line_callback:
            line_callback(xpixel_index + 1)
    return im

def get_heatmap_image(M, entry_to_rgb):
    """
    The input is meant to be an entrywise squared correlation matrix.
    The indices are expected to have been reordered already.
    @param M: a matrix with entries between zero and one
    @param entry_to_rgb: a function that returns an RGB color given an entry of M
    @return: a PIL image
    """
    n = len(M)
    im = Image.new('RGB', (n, n), 'red')
    for i in range(n):
        for j in range(n):
            rgb = entry_to_rgb(M[i][j])
            im.putpixel((i, j), rgb)
    return im

def get_block_heatmap_image(M, entry_to_rgb, blockwidth=3):
    """
    Use blocks larger than pixels.
    The input is meant to be an entrywise squared correlation matrix.
    The indices are expected to have been reordered already.
    @param M: a matrix with entries between zero and one
    @param blockwidth: the width of one side of a block; defaults to a single pixel
    @param entry_to_rgb: a function that returns an RGB color given an entry of M
    @return: a PIL image
    """
    #TODO use real PIL drawing functions to draw the blocks
    n = len(M)
    im = Image.new('RGB', (n*blockwidth, n*blockwidth), 'red')
    for i in range(n):
        for j in range(n):
            rgb = entry_to_rgb(M[i][j])
            for u in range(blockwidth):
                for v in range(blockwidth):
                    x = i*blockwidth + u
                    y = j*blockwidth + v
                    im.putpixel((x, y), rgb)
    return im

def get_heatmap(RoR, filename):
    """
    The input is meant to be an entrywise squared correlation matrix.
    The indices are expected to have been reordered already.
    @param RoR: a matrix with entries between zero and one
    @param filename: where to put the image
    """
    #FIXME this function is pretty useless
    f = gradient.squared_correlation_to_rgb
    im = get_heatmap_image(RoR, f)
    fout = open(filename, 'wb')
    im.save(fout)
    fout.close()

class DendrogramImager:
    """
    This is a vehicle for the dendrogram line drawing function.
    A PIL RGB image is created.
    This class is to bitmaps as dendro.AsciiArt is to ascii art.
    """

    def __init__(self, root):
        """
        Initialization creates the dendrogram image.
        @param: the root of a tree
        """
        breadth_gap = 2
        height_gap = 3
        image_height = dendro.get_dendrogram_breadth(root, breadth_gap)
        image_width = dendro.get_dendrogram_height(root, height_gap)
        self.im = Image.new('RGB', (image_width, image_height), 'white')
        dendro.draw_dendrogram(root, breadth_gap, height_gap, self.on_draw_line)

    def on_draw_line(self, line):
        """
        @param line: each endpoint is a (breadth_offset, height_offset) pair
        """
        #TODO draw a line using PIL instead of putting the pixels one by one
        ((a, b), (c, d)) = line
        if a == c:
            for height_offset in range(min(b, d), max(b, d) + 1):
                self.im.putpixel((height_offset, a), (0, 0, 0))
        elif b == d:
            for breadth_offset in range(min(a, c), max(a, c) + 1):
                self.im.putpixel((b, breadth_offset), (0, 0, 0))


def get_heatmap_with_dendrogram(M, root, entry_to_rgb, filename):
    """
    The input matrix is meant to be an entrywise squared correlation matrix.
    The indices of M should have been ordered conformantly with the dendrogram order.
    @param M: a matrix with entries between zero and one
    @param root: the root of the tree used to create the dendrogram
    @param entry_to_rgb: a function that returns an RGB color given an entry of M
    @param filename: where to put the image
    """
    #TODO add an option for the side where the dendrogram goes
    left = True
    # draw the dendrogram
    imager = DendrogramImager(root)
    dendrogram_image = imager.im
    # draw the heatmap
    heatmap_image = get_block_heatmap_image(M, entry_to_rgb)
    # paste the heatmap and the dendrogram into the same image
    n = len(M)
    image_width = n*3 + dendrogram_image.size[0] + 2
    image_height = n*3
    composite_image = Image.new('RGB', (image_width, image_height), 'white')
    composite_image.paste(heatmap_image, (0, 0))
    composite_image.paste(dendrogram_image, (n*3 + 1, 1))
    fout = open(filename, 'wb')
    composite_image.save(fout)
    fout.close()

def get_newick_order(s):
    """
    Extract a permutation from a newick string.
    @param s: a newick string
    @return: a permutation of the first N nonnegative integers
    """
    # first change every nondigit value in s to a space
    lex = []
    for c in s:
        if c not in '0123456789':
            lex.append(' ')
        else:
            lex.append(c)
    spaced_s = ''.join(lex)
    return [int(x) for x in spaced_s.split()]


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

    def test_dendrogram_imager(self):
        filename = 'dendrogram-test.png'
        root = mtree.create_tree([[0, [1, 2]], 3, [4, 5, 6]])
        imager = DendrogramImager(root)
        fout = open(filename, 'wb')
        imager.im.save(fout)
        fout.close()

    def test_heatmap_with_dendrogram(self):
        root = mtree.create_tree([[0, [1, 2]], 3, [4, 5, 6]])
        M = np.random.random((7, 7))
        R = np.corrcoef(M)
        # draw the correlation heatmap
        filename = 'r-test.png'
        f = gradient.correlation_to_rgb
        get_heatmap_with_dendrogram(R, root, f, filename)
        # draw the squared correlation heatmap
        filename = 'rr-test.png'
        RoR = R*R
        f = gradient.squared_correlation_to_rgb
        get_heatmap_with_dendrogram(RoR, root, f, filename)


if __name__ == '__main__':
    unittest.main()
