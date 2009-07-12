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


class TreeData:
    """
    This is used to maintain state whose scope is the construction of a tree.
    Please do not move this stuff into the mtree itself.
    """

    def __init__(self):
        """
        Initialize stuff that is meaningful per tree construction.
        """
        # leaves representing outgroups have negative labels
        self._next_outgroup_label = -1
        # maintain a map from labels to nodes
        self.label_to_node = {}

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


def build_tree(U_in, S_in, ordered_labels, tree_data):
    """
    Get the root of an mtree reconstructed from the transformed data.
    The input matrix U will be freed (deleted) by this function.
    @param U_in: part of the laplacian sqrt obtained by svd
    @param S_in: another part of the laplacian sqrt obtained by svd
    @param ordered_labels: a list of labels conformant with rows of U
    @param tree_data: state whose scope is the construction of the tree
    @return: an mtree rooted at a degree 2 vertex unless the input matrix has 3 rows
    """
    p, n = U_in.shape
    if p < 3:
        raise ValueError('expected the input matrix to have at least three rows')
    # look for an informative split
    index_split = None
    if p > 3:
        # the signs of v match the signs of the fiedler vector
        v = khorr.get_fiedler_vector(U_in, S_in)
        index_split = splitbuilder.eigenvector_to_split(v)
        # if the split is degenerate then don't use it
        if min(len(x) for x in index_split) < 2:
            index_split = None
    # if no informative split was found then create a degenerate tree
    if not index_split:
        root = mtree.create_tree(ordered_labels)
        for node in root.preorder():
            if node.has_label():
                tree_data.add_node(node)
        return root
    # get the indices defined by the split
    a, b = tuple(list(sorted(x)) for x in index_split)
    # Create two new matrices.
    # Be somewhat careful to not create lots of intermediate matrices
    A = np.zeros((len(a)+1, n))
    B = np.zeros((len(b)+1, n))
    for i, index in enumerate(a):
        A[i] = U_in[index] * S_in
    for i, index in enumerate(b):
        B[i] = U_in[index] * S_in
    A_outgroup = np.sum(B, 0)
    B_outgroup = np.sum(A, 0)
    A[-1] = A_outgroup
    B[-1] = B_outgroup
    # delete the old matrix
    del U_in
    # recursively construct the subtrees
    subtrees = []
    for selection, complement, summed_L_sqrt in ((a,b,A),(b,a,B)):
        # record the outgroup label for this subtree
        outgroup_label = tree_data.decrement_outgroup_label()
        # create the ordered list of labels corresponding to leaves of the subtree
        next_ordered_labels = [ordered_labels[i] for i in selection]
        next_ordered_labels.append(outgroup_label)
        # get the criterion matrix for the next iteration
        U, S, VT = np.linalg.svd(summed_L_sqrt, full_matrices=0)
        # delete matrices that are no longer useful
        del summed_L_sqrt
        # build the tree recursively
        root = build_tree(U, S, next_ordered_labels, tree_data)
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

    def test_kh_dataset(self):
        X = splitbuilder.get_data('kh-dataset.txt')
        U, S = khorr.data_to_laplacian_sqrt(X)
        tree_data = TreeData()
        root = build_tree(U, S, range(len(U)), tree_data)
        #print
        #print root.get_newick_string()
        #print

    def test_fivetimes_dataset(self):
        X = splitbuilder.get_data('fivetimes2.txt')
        X = X.T
        U, S = khorr.data_to_reduced_laplacian_sqrt(X)
        tree_data = TreeData()
        root = build_tree(U, S, range(len(U)), tree_data)
        for node in root.preorder():
            if node.has_label():
                node.label += 1
        print
        print root.get_newick_string()
        print


if __name__ == '__main__':
    unittest.main()
