"""
This module is for drawing a dendrogram.

Branch lengths are disregarded.
This could be made more object oriented.
"""

import unittest

import mtree

def get_dendrogram_breadth(root, breadth_gap):
    """
    @param root: the root of a tree
    @param breadth_gap: greater to draw broad trees
    @return: a number of layout units (e.g. pixels, characters)
    """
    return (len(root.ordered_tips()) - 1) * (breadth_gap + 1) + 1

def get_dendrogram_height(root, height_gap):
    """
    @param root: the root of a tree
    @param height_gap: greater to draw tall trees
    @return: a number of layout units (e.g. pixels, characters)
    """
    return root.height() * (height_gap + 1) + 1

def draw_dendrogram(root, breadth_gap, height_gap, draw_line):
    """
    This is the main interface for drawing the dendrogram.
    The line drawing function takes a pair of points.
    Each point is a (breadth_offset, height_offset) pair.
    @param root: the root of a tree
    @param breadth_gap: greater to draw broad trees
    @param height_gap: greater to draw tall trees
    @param draw_line: a function that is called to draw a line
    """
    # compute heights
    id_to_height = build_heights(root, {})
    # compute the breadth offsets of the tips of the tree
    id_to_breadth_offset = {}
    for i, tip in enumerate(root.ordered_tips()):
        breadth_offset = i * (breadth_gap + 1)
        id_to_breadth_offset[id(tip)] = breadth_offset
    # recursively compute breadth offsets for non-tip nodes
    build_breadth_offsets(root, id_to_breadth_offset)
    # compute height offsets
    id_to_height_offset = {}
    for node in root.preorder():
        height_offset = id_to_height[id(node)] * (height_gap + 1)
        id_to_height_offset[id(node)] = height_offset
    # draw the broad lines
    for node in root.preorder():
        if node.children:
            height_offset = id_to_height_offset[id(node)]
            breadth_offsets = [id_to_breadth_offset[id(p)] for p in node.children]
            a = (min(breadth_offsets), height_offset)
            b = (max(breadth_offsets), height_offset)
            draw_line((a, b))
    # draw the tall lines
    for node in root.preorder():
        if node.parent:
            breadth_offset = id_to_breadth_offset[id(node)]
            height_offsets = [id_to_height_offset[id(p)] for p in (node, node.parent)]
            a = (breadth_offset, min(height_offsets))
            b = (breadth_offset, max(height_offsets))
            draw_line((a, b))

def draw_tall_dendrogram(root, breadth_gap, height_gap, draw_line):
    """
    This function makes ugly dendrograms.
    The line drawing function takes a pair of points.
    Each point is a (breadth_offset, height_offset) pair.
    @param root: the root of a tree
    @param breadth_gap: greater to draw broad trees
    @param height_gap: greater to draw tall trees
    @param draw_line: a function that is called to draw a line
    """
    # compute distances to the root
    id_to_topo_depth = build_topo_depths(root, 0, {})
    max_topo_depth = root.height()
    # compute the breadth offsets of the tips of the tree
    id_to_breadth_offset = {}
    for i, tip in enumerate(root.ordered_tips()):
        breadth_offset = i * (breadth_gap + 1)
        id_to_breadth_offset[id(tip)] = breadth_offset
    # recursively compute breadth offsets for non-tip nodes
    build_breadth_offsets(root, id_to_breadth_offset)
    # compute height offsets
    id_to_height_offset = {}
    for node in root.preorder():
        if node.children:
            topo_depth = id_to_topo_depth[id(node)]
            height_offset = (max_topo_depth - topo_depth) * (height_gap + 1)
        else:
            height_offset = 0
        id_to_height_offset[id(node)] = height_offset
    # draw the broad lines
    for node in root.preorder():
        if node.children:
            height_offset = id_to_height_offset[id(node)]
            breadth_offsets = [id_to_breadth_offset[id(p)] for p in node.children]
            a = (min(breadth_offsets), height_offset)
            b = (max(breadth_offsets), height_offset)
            draw_line((a, b))
    # draw the tall lines
    for node in root.preorder():
        if node.parent:
            breadth_offset = id_to_breadth_offset[id(node)]
            height_offsets = [id_to_height_offset[id(p)] for p in (node, node.parent)]
            a = (breadth_offset, min(height_offsets))
            b = (breadth_offset, max(height_offsets))
            draw_line((a, b))

def build_heights(node, id_to_height):
    """
    This recursive function build distances to the furthest leaf.
    @param node: a tree node
    @param id_to_height: the dictionary of node ids to heights
    """
    for child in node.children:
        build_heights(child, id_to_height)
    height = 1 + max([-1] + [id_to_height[id(c)] for c in node.children])
    id_to_height[id(node)] = height
    return id_to_height

def build_topo_depths(node, topo_depth, id_to_topo_depth):
    """
    This recursive function builds distances to the root.
    @param node: a tree node
    @param topo_depth: the distance from the node to the root
    @param id_to_topo_depth: a map from tree node id to distance from root
    """
    id_to_topo_depth[id(node)] = topo_depth
    for child in node.children:
        build_topo_depths(child, topo_depth+1, id_to_topo_depth)
    return id_to_topo_depth

def build_breadth_offsets(node, id_to_breadth_offset):
    """
    This recursive function builds breadth offsets.
    Tips of the tree are expected to have known breadth offsets.
    @param node: a tree node
    @param id_to_breadth_offset: a map from tree node id to breadth offset
    @return: the breadth offset
    """
    if not node.children:
        return id_to_breadth_offset[id(node)]
    child_offsets = [build_breadth_offsets(child, id_to_breadth_offset) for child in node.children]
    offset = (max(child_offsets) + min(child_offsets)) / 2
    id_to_breadth_offset[id(node)] = offset
    return offset


class AsciiArt:
    """
    This is a vehicle for the drawline callback.
    """

    def __init__(self, root, draw_function=draw_dendrogram):
        self.root = root
        self.draw_function = draw_function

    def on_draw_line(self, line):
        """
        @param line: each endpoint is a (breadth_offset, height_offset) pair
        """
        ((a, b), (c, d)) = line
        if a == c:
            for height_offset in range(min(b, d), max(b, d) + 1):
                self.art[a][height_offset] = '-'
        elif b == d:
            for breadth_offset in range(min(a, c), max(a, c) + 1):
                self.art[breadth_offset][b] = '|'
        for breadth_offset, height_offset in line:
            self.art[breadth_offset][height_offset] = '+'

    def __str__(self):
        breadth_gap = 2
        height_gap = 3
        breadth = get_dendrogram_breadth(self.root, breadth_gap)
        height = get_dendrogram_height(self.root, height_gap)
        self.art = [['.']*height for i in range(breadth)]
        self.draw_function(self.root, breadth_gap, height_gap, self.on_draw_line)
        return '\n'.join(''.join(row) for row in self.art)


g_expected_tall_art = """+-------+....
........|....
........+---+
+---+...|...|
....+---+...|
....|.......|
+---+.......|
............|
............|
+-----------+
............|
............|
+-------+...|
........|...|
........|...|
+-------+---+
........|....
........|....
+-------+...."""

g_expected_short_art = """+-------+....
........|....
........+---+
+---+...|...|
....+---+...|
....|.......|
+---+.......|
............|
............|
+-----------+
............|
............|
+---+.......|
....|.......|
....|.......|
+---+-------+
....|........
....|........
+---+........"""


class TestMe(unittest.TestCase):

    def test_dendrogram(self):
        root = mtree.create_tree([[0, [1, 2]], 3, [4, 5, 6]])
        # make a tall tree
        observed_tall_art = str(AsciiArt(root, draw_tall_dendrogram))
        self.assertEqual(observed_tall_art, g_expected_tall_art)
        # make a short tree
        observed_short_art = str(AsciiArt(root))
        self.assertEqual(observed_short_art, g_expected_short_art)

if __name__ == '__main__':
    unittest.main()


