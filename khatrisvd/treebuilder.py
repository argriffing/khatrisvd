"""
This module is supposed to actually build the tree from the data.

Instead of repeatedly splitting a distance matrix,
the splitting function will repeatedly split a rectangular matrix.
"""

import logging
import unittest

import numpy as np

import mtree
import khorr
import splitbuilder
import util

def build_naive_divisive_tree_with_early_termination(X):
    """
    Build a tree by recursively splitting gene sets without regard to outgrouping.
    @param X: a data matrix, preferably with more rows than columns
    @return: the root of a tree
    """
    return build_naive_divisive_tree(X, early_termination=True)

def build_naive_divisive_tree(X, early_termination=False):
    """
    Build a tree by recursively splitting gene sets without regard to outgrouping.
    This method is less naive than building from a single split,
    and is more naive than splitting with outgrouping.
    @param X: a data matrix, preferably with more rows than columns
    @param early_termination: True iff clustering stops when a split is degenerate
    @return: the root of a tree
    """
    p, n = X.shape
    Z = khorr.get_standardized_matrix(X)
    Z = khorr.standardized_to_augmented_C(Z)
    boxed_Z = [Z]
    del Z
    return _build_naive_divisive_tree_helper(boxed_Z, range(p), early_termination)

def _build_naive_divisive_tree_helper(boxed_Z, ordered_labels, early_termination=False):
    """
    Try to be somewhat memory efficient because Z can be huge.
    @param boxed_Z: a standardized data matrix, boxed so it can be deleted
    @param ordered_labels: integer labels conformant to rows of Z
    @param early_termination: True iff clustering stops when a split is degenerate
    @return: the root of a tree
    """
    if len(boxed_Z) != 1:
        raise ValueError('expected the input matrix to be boxed for deletion')
    Z = boxed_Z[0]
    if len(Z) != len(ordered_labels):
        raise ValueError('the input labels are incompatible with the input matrix')
    p = len(ordered_labels)
    # define the root
    root = mtree.Node()
    # deal with a degenerate split
    if p == 1:
        root.label = ordered_labels[0]
        return root
    # get the eigenvector whose loadings will be used to split the matrix
    Z = util.get_column_centered_matrix(Z)
    U, S, VT = np.linalg.svd(Z, full_matrices=0)
    v = khorr.get_dominant_vector(U, S)
    del U
    del VT
    # split the matrix
    stack = []
    index_split = splitbuilder.eigenvector_to_split(v)
    # if we are doing early termination and the split is degenerate then we are done
    if early_termination and min(len(x) for x in index_split) < 2:
        for loading, row_index in sorted((x, i) for i, x in enumerate(v)):
            child = mtree.Node()
            child.label = ordered_labels[row_index]
            root.add_child(child)
        return root
    for selection_set in index_split:
        selection = list(sorted(selection_set))
        # define the next standardized (but not column centered) matrix
        next_matrix = np.vstack(row for i, row in enumerate(Z) if i in selection_set)
        # define the next ordered labels
        next_ordered_labels = [ordered_labels[i] for i in selection]
        # add to the stack
        stack.append([next_matrix, next_ordered_labels])
    # we no longer need the Z matrix
    del boxed_Z[0]
    del Z
    # build the tree
    while stack:
        next_matrix, next_ordered_labels = stack.pop()
        next_boxed_Z = [next_matrix]
        del next_matrix
        child = _build_naive_divisive_tree_helper(next_boxed_Z, next_ordered_labels, early_termination)
        root.add_child(child)
    return root

def build_single_split_correlation_tree(X):
    """
    This is a naive method to use as a control.
    I think it is equivalent to eigengene computations.
    @param X: a data matrix, preferably with more rows than columns
    @return: the root of a tree
    """
    return build_single_split_tree(X, use_squared_correlation=False)

def build_single_split_tree(X, use_squared_correlation=True):
    """
    Get the root of an mtree reconstructed from the transformed data.
    Note that only the dominant singular vector is required.
    This may be faster to get than the entire SVD.
    This method is naive compared to build_tree.
    With the use_squared_correlation options disabled, it is even more naive.
    @param X: a data matrix, preferably with more rows than columns
    @param use_squared_correlation: True for squared correlation, False for correlation
    @return: the root of a tree
    """
    # get the eigenvector whose loadings will be used to split and order the rows
    logging.debug('creating the standardized matrix')
    Z = khorr.get_standardized_matrix(X)
    if use_squared_correlation:
        logging.debug('creating the augmented matrix')
        Z = khorr.standardized_to_augmented_C(Z)
    logging.debug('creating the column centered matrix')
    W = util.get_column_centered_matrix(Z)
    logging.debug('manually cleaning up old matrices')
    del Z
    logging.debug('doing a singular value decomposition')
    U, S, VT = np.linalg.svd(W, full_matrices=0)
    logging.debug('getting the dominant eigenvector')
    v = khorr.get_dominant_vector(U, S)
    # account for values near zero, using the same criterion as in splitbuilder
    epsilon = 1e-14
    vprime = [0.0 if abs(x) < epsilon else x for x in v]
    # start making a tree from the eigenvector
    root = mtree.Node()
    neg_child = mtree.Node()
    pos_child = mtree.Node()
    root.add_child(neg_child)
    root.add_child(pos_child)
    for loading, row_index in sorted((x, i) for i, x in enumerate(vprime)):
        grandchild = mtree.Node()
        grandchild.label = row_index
        if loading > 0:
            pos_child.add_child(grandchild)
        else:
            neg_child.add_child(grandchild)
    return root

def build_tree(X):
    """
    Get the root of an mtree reconstructed from the transformed data.
    @param X: a data matrix, preferably with more rows than columns
    @return: the root of a tree
    """
    U, S = khorr.data_to_reduced_laplacian_sqrt(X)
    tree_data = TreeData()
    # The prancing around with U is so that it can be deleted in the called function.
    boxed_U = [U]
    del U
    root = build_tree_helper(boxed_U, S, range(len(X)), tree_data)
    sort_tree(root)
    return root


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


def build_tree_helper(boxed_U_in, S_in, ordered_labels, tree_data):
    """
    Get the root of an mtree reconstructed from the transformed data.
    The input matrix U will be freed (deleted) by this function.
    @param boxed_U_in: part of the laplacian sqrt obtained by svd
    @param S_in: another part of the laplacian sqrt obtained by svd
    @param ordered_labels: a list of labels conformant with rows of U
    @param tree_data: state whose scope is the construction of the tree
    @return: an mtree rooted at a degree 2 vertex unless the input matrix has 3 rows
    """
    # take U_in out of the box
    if len(boxed_U_in) != 1:
        raise ValueError('expected a 2d array as the only element of a list')
    U_in = boxed_U_in[0]
    shape = U_in.shape
    if len(shape) != 2:
        raise valueError('expected a 2d array as the only element of a list')
    p, n = shape
    if p < 3 or n < 3:
        raise ValueError('expected the input matrix to have at least three rows and columns')
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
    # delete the two references to the old matrix
    del U_in
    del boxed_U_in[0]
    # recursively construct the subtrees
    subtrees = []
    stack = [[b,a,B], [a,b,A]]
    # delete non-stack references to partial matrices
    del A
    del B
    # process the partial matrices
    while stack:
        selection, complement, summed_L_sqrt = stack.pop()
        # record the outgroup label for this subtree
        outgroup_label = tree_data.decrement_outgroup_label()
        # create the ordered list of labels corresponding to leaves of the subtree
        next_ordered_labels = [ordered_labels[i] for i in selection]
        next_ordered_labels.append(outgroup_label)
        # get the criterion matrix for the next iteration
        U, S, VT = np.linalg.svd(summed_L_sqrt, full_matrices=0)
        del VT
        # delete matrices that are no longer useful
        del summed_L_sqrt
        # build the tree recursively
        boxed_U = [U]
        del U
        root = build_tree_helper(boxed_U, S, next_ordered_labels, tree_data)
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

def sort_tree(root):
    """
    Enforce an ordering on the child trees.
    The children of a node will be sorted according to their size,
    where their size is defined by the number of leaves in their subtree.
    @param root: the root of an mtree
    @return: the number of leaves in the subtree defined by the root.
    """
    if not root.children:
        return 1
    count_child_pairs = list(sorted((sort_tree(child), child) for child in root.children))
    # reorder the children
    root.children = [child for count, child in count_child_pairs]
    # return the number of leaves
    return sum(count for count, child in count_child_pairs)


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
        X = util.file_to_whitespace_separated_matrix('kh-dataset.txt')
        root = build_tree(X)
        #print
        #print root.get_newick_string()
        #print

    def test_fivetimes_dataset(self):
        X = util.file_to_whitespace_separated_matrix('fivetimes2.txt')
        # this file was generated transposed; oops
        X = X.T
        # build the tree
        print X.shape
        root = build_tree(X)
        for node in root.preorder():
            if node.has_label():
                node.label += 1
        print
        print root.get_newick_string()
        print


if __name__ == '__main__':
    unittest.main()
