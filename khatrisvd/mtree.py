"""
This module will deal with lightweight multifurcating trees.

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

    def has_label(self):
        return self.label is not None

    def preorder(self):
        """
        @return: a list of vertices in a partial order
        """
        vertices = [self]
        i = 0
        while i < len(vertices):
            vertices.extend(vertices[i].children)
            i += 1
        return vertices

    def remove_child(self, child):
        """
        This function breaks the tree into two trees.
        @param child: a child that will become the root of its own tree
        """
        child.parent = None
        self.children.remove(child)

    def add_child(self, node):
        """
        This function assumes that the node is not already a child.
        @param node: the node to be added as a child
        """
        node.parent = self
        self.children.append(node)

    def is_root(self):
        """
        @return: True iff this node is a root
        """
        return not self.parent

    def gen_path_to_root(self):
        """
        Yield nodes on the way up to the root.
        """
        cur = self
        while cur:
            yield cur
            cur = cur.parent

    def get_path_from_root(self):
        """
        Return the list of nodes on the path from the root to this node.
        """
        p = list(self.gen_path_to_root())
        p.reverse()
        return p

    def reroot(self):
        """
        Reroot the tree to this node.
        """
        if self.is_root():
            return
        p = self.get_path_from_root()
        for root, child in zip(p[:-1], p[1:]):
            root.remove_child(child)
            child.add_child(root)

    def remove(self):
        """
        Remove the node of degree 1 or 2.
        @return: the formerly neighboring node that is now closest to the root
        """
        if self.degree() not in (1, 2):
            raise ValueError('only nodes of degree 1 or 2 can be removed')
        if not self.parent:
            self.children[0].reroot()
        parent = self.parent
        if self.children:
            parent.add_child(self.children[0])
        parent.remove_child(self)
        self.children = []
        return parent

    def degree(self):
        """
        @return: the degree of the node
        """
        d = len(self.children)
        if self.parent:
            d += 1
        return d

    def get_newick_string(self):
        """
        @return: a newick string
        """
        return self.get_newick_substring() + ';'

    def get_newick_substring(self):
        """
        @return: part of a newick string corresponding to the subtree
        """
        name = str(self.label) if self.has_label() else ''
        if not self.children:
            return name
        else:
            return '(' + ', '.join(child.get_newick_substring() for child in self.children) + ')' + name


def create_tree(term):
    """
    This recursive function builds a tree from builtin types.
    @param term: a label or a list of lists and labels
    @return: the root of the new tree
    """
    root = Node()
    try:
        for subterm in term:
            root.add_child(create_tree(subterm))
    except TypeError, e:
        root.label = term
    return root


class TestMe(unittest.TestCase):

    def test_creation(self):
        root = create_tree([[0, 1, 2], 3, [4, 5, 6]])

    def test_newick(self):
        root = create_tree([[0, [1, 2]], 3, [4, 5, 6]])
        observed = root.get_newick_string()
        expected = '((0, (1, 2)), 3, (4, 5, 6));'
        self.assertEquals(observed, expected)

    def test_preorder(self):
        root = create_tree([[[0, 1]]])
        nodes = root.preorder()
        self.assertEqual(len(nodes), 5)
        expected = [False, False, False, True, True]
        observed = [v.has_label() for v in nodes]
        self.assertEqual(expected, observed)
        expected = [1, 2, 3, 1, 1]
        observed = [v.degree() for v in nodes]
        self.assertEqual(expected, observed)

    def test_reroot(self):
        root = create_tree([[[0, 1]]])
        root.reroot()
        nodes = root.preorder()
        nodes[-1].reroot()
        observed = nodes[-1].get_newick_string()
        expected = '((0, ()))1;'
        self.assertEqual(expected, observed)
        root.reroot()
        observed = root.get_newick_string()
        expected = '(((0, 1)));'
        self.assertEqual(expected, observed)

    def test_remove(self):
        # remove the root of a degenerate tree
        root = create_tree([[[0, 1]]])
        root = root.remove()
        expected = '((0, 1));'
        observed = root.get_newick_string()
        self.assertEqual(expected, observed)
        # remove a node after the root of a degenerate tree
        root = create_tree([[[0, 1]]])
        root = root.children[0].remove()
        expected = '((0, 1));'
        observed = root.get_newick_string()
        self.assertEqual(expected, observed)
        # remove the root of a bifurcating tree
        root = create_tree([[0,1],[2,3]])
        root = root.remove()
        expected = '(0, 1, (2, 3));'
        observed = root.get_newick_string()
        self.assertEqual(expected, observed)
        # trying to remove a degree three node fails
        self.assertRaises(ValueError, root.remove)


if __name__ == '__main__':
    unittest.main()
