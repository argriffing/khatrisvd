"""
This module is supposed to actually build the tree from the data.

Instead of repeatedly splitting a distance matrix,
the splitting function will repeatedly split a rectangular matrix.
"""


import mtree
import khorr
import splitbuilder

import numpy as np

import unittest


def update_svd(L_sqrt, row_index_complement):
    """
    Update the matrix sqrt of Laplacian by summing rows of the removed indices.
    The summed row goes at the end of the output matrix.
    @param L_sqrt: the matrix square root of a Laplacian
    @param row_index_complement: the set of row indices in the outgroup
    @return: an updated matrix sqrt of Laplacian
    """
    n = len(L_sqrt)
    row_index_selection = set(range(n)) - row_index_complement
    ordered_ingroup_rows = [L_sqrt[i] for i in sorted(row_index_selection)]
    ordered_outgroup_rows = [L_sqrt[i] for i in sorted(row_index_complement)]
    return np.vstack(ordered_ingroup_rows + [sum(ordered_outgroup_rows)])


class TreeData:
    """
    This is used to maintain state whose scope is the construction of a tree.
    Please do not move this stuff into the mtree itself.
    """

    def __init__(self, split_function, update_function):
        """
        Initialize stuff that is meaningful per tree construction.
        @param split_function: takes a matrix sqrt of Laplacian and returns an index split
        @param update_function: takes a matrix sqrt of Laplacian and an index subset and returns a matrix sqrt of Laplacian
        """
        # leaves representing outgroups have negative labels
        self._next_outgroup_label = -1
        # maintain a map from labels to nodes
        self.label_to_node = {}
        # save the split and update functions
        self.split_function = split_function
        self.update_function = update_function

    def decrement_outgroup_label(self):
        """
        Return a tree-scoped outgroup label with side effect.
        The side effect is to decrement the label that is on deck.
        @return: a label that is a negative integer
        """
        label = self._next_outgroup_label
        self._next_outgroup_label -= 1
        return label

    def add_node(self, node):
        """
        @param node: a node to be added to the label_to_node map
        """
        label = node.label
        if label in self.label_to_node:
            raise ValueError('the label of this node is already in the dictionary')
        self.label_to_node[label] = node

    def remove_node(self, node):
        """
        @param node: a node to be removed from the label_to_node map
        """
        del self.label_to_node[node.label]


def build_tree(L_sqrt, ordered_labels, tree_data):
    """
    Get the root of an mtree reconstructed from the transformed data.
    @param L_sqrt: the matrix square root of a Laplacian
    @param ordered_labels: a list of labels conformant with rows of L_sqrt
    @param tree_data: state whose scope is the construction of the tree
    @return: an mtree rooted at a degree 2 vertex unless the input matrix has 3 rows
    """
    #TODO make this function iterative instead of recursive
    n = len(L_sqrt)
    if n < 3:
        raise ValueError('expected the input matrix to have at least three rows')
    # look for an informative split
    index_split = None
    if n > 3:
        try:
            index_split = tree_data.split_function(L_sqrt)
        except splitbuilder.DegenerateSplitException, e:
            pass
    # if no informative split was found then create a degenerate tree
    if not index_split:
        root = mtree.create_tree(ordered_labels)
        for node in root.preorder():
            if node.has_label():
                tree_data.add_node(node)
        return root
    # recursively construct the subtrees
    subtrees = []
    a, b = index_split
    for index_selection, index_complement in ((a,b),(b,a)):
        # record the outgroup label for this subtree
        outgroup_label = tree_data.decrement_outgroup_label()
        # create the ordered list of labels corresponding to leaves of the subtree
        next_ordered_labels = [ordered_labels[i] for i in sorted(index_selection)]
        next_ordered_labels.append(outgroup_label)
        # create the next matrix with rows conformant to the ordered labels
        next_L_sqrt = tree_data.update_function(L_sqrt, index_complement)
        # recursively construct the subtree
        root = build_tree(next_L_sqrt, next_ordered_labels, tree_data)
        # if the root is degree 2 then remove the root node
        if root.degree() == 2:
            root = root.remove()
        # root the tree at the outgroup node
        root = tree_data.label_to_node[outgroup_label]
        root.reroot()
        # we don't need the outgroup label anymore
        tree_data.remove_node(root)
        # we can also remove the label from the outgroup node itself
        root.label = None
        # save the properly rooted subtree
        subtrees.append(root)
    # connect the two subtrees at their roots
    left_root, right_root = subtrees
    right_root = right_root.remove()
    left_root.add_child(right_root)
    return left_root


class TestMe(unittest.TestCase):

    def assertAllClose(self, A, B, msg=None):
        """
        Assert that two matrices are almost entrywise identical.
        @param A: a matrix
        @param B: a matrix
        @param msg: shown when the assertion fails
        """
        self.assertTrue(np.allclose(A, B), msg=msg)

    def test_update_svd(self):
        index_set = set([1, 2])
        M = np.array([[1,2,3,4],[5,6,7,8],[9,10,11,12],[13,14,15,16]])
        expected = np.array([[1,2,3,4],[13,14,15,16],[14,16,18,20]])
        observed = update_svd(M, index_set)
        self.assertAllClose(observed, expected)

    def test_kh_dataset(self):
        X = splitbuilder.get_data('kh-dataset.txt')
        L_sqrt = khorr.data_to_laplacian_sqrt(X)
        tree_data = TreeData(splitbuilder.split_svd, update_svd)
        root = build_tree(L_sqrt, range(len(X)), tree_data)
        #print
        #print root.get_newick_string()
        #print

    def test_fivetimes_dataset(self):
        X = splitbuilder.get_data('fivetimes2.txt')
        X = X.T
        L_sqrt = khorr.data_to_laplacian_sqrt(X)
        tree_data = TreeData(splitbuilder.split_svd, update_svd)
        root = build_tree(L_sqrt, range(len(X)), tree_data)
        for node in root.preorder():
            if node.has_label():
                node.label += 1
        print
        print root.get_newick_string()
        print


if __name__ == '__main__':
    unittest.main()
