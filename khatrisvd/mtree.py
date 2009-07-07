"""
This module will deal with multifurcating trees.

The tree will not have branch lengths.
Each tip of the tree will be labeled with an integer.
The root node of the tree will be used instead
of a separate enclosing tree class.
Trees are expected to get huge enough that
you probably don't want to do anything N^2.
Also iterative algorithms should probably be used instead
of recursive algorithms because python
is not so great for recursion.
"""

import unittest


class Node:

    def __init__(self):
        self.label = None
        self.parent = None
        self.children = []

    def reroot(self):
        """
        Reroot the tree at this node.
        """
        #TODO do i even need this
        pass

    def foo(self):
        pass


#TODO add tests


if __name__ == '__main__':
    unittest.main()
