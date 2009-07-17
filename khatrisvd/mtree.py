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

    def height(self):
        """
        @return: the height of the tree
        """
        if not self.children:
            return 0
        return max(child.height() for child in self.children) + 1

    def ordered_tips(self):
        """
        @return: a list of vertices in an order suitable for a dendrogram
        """
        if not self.children:
            return [self]
        else:
            arr = []
            for child in self.children:
                arr.extend(child.ordered_tips())
            return arr

    def ordered_labels(self):
        """
        @return: a list of labels in an order suitable for a dendrogram
        """
        return [node.label for node in self.ordered_tips()]

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

    def remove_stable_degree_two_non_root(self):
        """
        Remove a useless node.
        @return: the node that was the parent
        """
        if not self.parent:
            raise ValueError('this function does not work on the root')
        if len(self.children) != 1:
            raise ValueError('this function does not work with more than one child')
        parent = self.parent
        child = self.children[0]
        self.children = []
        self.parent = None
        parent.children = [(child if c is self else c) for c in parent.children]
        child.parent = parent
        return parent

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

def _newick_to_tokens(newick):
    """
    The newick string must be absurdly simple.
    @return: a list of tokens
    """
    tokens = []
    int_string = ''
    digits = '0123456789'
    punctuation = '(),;'
    for c in newick:
        if c in digits:
            int_string += c
        else:
            if int_string:
                tokens.append(int(int_string))
                int_string = ''
            if c in punctuation:
                tokens.append(c)
    return tokens

def _tokens_to_subtree(tokens, initial_offset):
    """
    The first indexed token should be an open parenthesis.
    The last indexed token should be after a close parenthesis.
    @param tokens: the list of tokens
    @param offset: the offset into the list of tokens
    @return: an offset and a subtree
    """
    root = Node()
    offset = initial_offset
    if type(tokens[offset]) is int:
        # handle the degenerate case
        root.label = tokens[offset]
    elif tokens[initial_offset] == '(':
        # skip past the open parenthesis
        offset += 1
        # go until a close parenthesis
        while tokens[offset] != ')':
            if tokens[offset] == '(' or type(tokens[offset]) is int:
                offset, child = _tokens_to_subtree(tokens, offset)
                root.add_child(child)
            else:
                offset += 1
    return offset + 1, root

def newick_file_to_mtree(newick_filename):
    fin = open(newick_filename)
    newick_string = fin.read()
    fin.close()
    return newick_to_mtree(newick_string)

def newick_to_mtree(newick):
    """
    The newick string must be absurdly simple.
    """
    tokens = _newick_to_tokens(newick)
    if tokens[-1] != ';':
        raise ValueError('the newick string must end with a semicolon')
    offset, root = _tokens_to_subtree(tokens, 0)
    return root

def _leaves_to_subtree_helper(root, required_ids):
    """
    @param root: the root of the template tree
    @param required_ids: clone these nodes in the template tree
    @return: a new node
    """
    if id(root) not in required_ids:
        raise ValueError('the root should be in the required id set')
    cloned_root = Node()
    if root.has_label():
        cloned_root.label = root.label
    for child in root.children:
        if id(child) in required_ids:
            cloned_root.add_child(_leaves_to_subtree_helper(child, required_ids))
    return cloned_root

def leaves_to_subtree(root, leaves):
    """
    None of the original nodes are in the returned subtree.
    @parm root: the root of a tree
    @param leaves: leaves of the tree
    @return: a copy of enough of the tree to reach the leaves
    """
    # get the set of required nodes
    required_ids = set()
    for node in leaves:
        current = node
        while current:
            if id(current) in required_ids:
                break
            required_ids.add(id(current))
            current = current.parent
    # clone the required subset of the tree
    cloned_root = _leaves_to_subtree_helper(root, required_ids)
    # remove non-root degree 2 nodes
    for node in cloned_root.preorder():
        if node.parent and node.degree() == 2:
            node.remove_stable_degree_two_non_root()
    # return the cloned root
    return cloned_root


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

    def test_ordered_labels(self):
        root = create_tree([[0, [1, 2]], 3, [4, 5, 6]])
        observed = root.ordered_labels()
        expected = range(7)
        self.assertEqual(expected, observed)

    def test_newick_to_tokens(self):
        newick = '(0, 1, (24, 3));'
        expected = ['(', 0, ',', 1, ',', '(', 24, ',', 3, ')', ')', ';']
        observed = _newick_to_tokens(newick)
        self.assertEqual(expected, observed)

    def test_newick_to_mtree(self):
        newick = '(0, 1, (24, 3));'
        root = newick_to_mtree(newick)
        observed = root.get_newick_string()
        expected = newick
        self.assertEqual(expected, observed)

    def test_leaves_to_subtree(self):
        root = create_tree([[0, [1, 2]], 3, [4, 5, 6]])
        leaves = [tip for tip in root.preorder() if tip.label in [1,2,3,4]]
        cloned = leaves_to_subtree(root, leaves)
        expected = '((1, 2), 3, 4);'
        observed = cloned.get_newick_string()
        self.assertEqual(expected, observed)


if __name__ == '__main__':
    unittest.main()
