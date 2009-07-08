"""
This module is supposed to actually build the tree from the data.

Instead of repeatedly splitting a distance matrix,
the splitting function will repeatedly split a rectangular matrix.
"""

#TODO some code has been copypasted from splitbuilder; remove duplicated code

import mtree
import khorr
import splitbuilder

import numpy as np

import unittest


class InvalidSpectralSplitException(Exception):
    """
    This exception is raised when a spectral split fails horribly.
    In particular, it is raised when an eigenvector based splitting method
    gets an eigenvector that does not define any split, even a degenerate split.
    """
    def __init__(self, L_sqrt):
        """
        @param L_sqrt: sqrt of Laplacian to save for posterity
        """
        Exception.__init__(self)
        self.L_sqrt = L_sqrt


class DegenerateSplitException(Exception):
    """
    This exception is raised when one taxon is split from the rest.
    This kind of split is not informative.
    """
    def __init__(self, index):
        """
        @param index: the degenerate index to blame later
        """
        Exception.__init__(self)
        self.index = index


def make_split(a, b):
    """
    This is a helper function.
    @param a: a sequence of hashable values
    @param b: a sequence of hashable values
    @return: a split
    """
    return frozenset((frozenset(a), frozenset(b)))

def set_to_string(my_set):
    """
    @param my_set: a sortable sequence
    @return: a string
    """
    return '{' + ', '.join(str(x) for x in sorted(my_set)) + '}'

def split_to_string(my_split):
    """
    @param my_split: a frozenset of some frozensets
    @return: a string
    """
    return set_to_string([set_to_string(x) for x in my_split])

def eigenvector_to_split(v, epsilon=1e-14):
    """
    This is a helper function.
    @param v: the signs of the loadings of this eigenvector define the split
    @param epsilon: negligible eigenvector loadings will be treated as zero
    @return: a split
    """
    vprime = [0.0 if abs(x) < epsilon else x for x in v]
    positive_indices = [i for i, x in enumerate(vprime) if x > 0]
    nonpositive_indices = [i for i, x in enumerate(vprime) if x <= 0]
    return make_split(positive_indices, nonpositive_indices)

def vmerge(label_sets, index_set):
    """
    @param label_sets: an ordered list of sets of disjoint integers
    @param index_set: a set of indices to merge
    @return: a stably merged list of label sets
    """
    # define the merged set of labels
    merged_labels = set()
    for i in index_set:
        merged_labels.update(label_sets[i])
    # get the min index which will be replaced
    min_index = min(index_set)
    # define the new list of label sets
    next_label_sets = []
    for i, label_set in enumerate(label_sets):
        if i == min_index:
            next_label_sets.append(merged_labels)
        elif i not in index_set:
            next_label_sets.append(label_set)
    return next_label_sets

def rmerge(M, index_set):
    """
    Merge rows indexed by an index set.
    The row indexed by the minimum index will be replaced by the sum of the rows.
    The other indexed rows will be removed.
    The order of rows indexed by indices not in the set will be stable.
    @param M: a matrix
    @param index_set: a set of row indices to merge
    """
    n = len(M)
    index_sets = vmerge([set([i]) for i in range(n)], index_set)
    n_small = len(index_sets)
    M_small = np.zeros((n, n_small))
    rows = []
    for i, set_i in enumerate(index_sets):
        rows.append(sum(M[k] for k in set_i))
    return np.array(rows)

def index_split_to_label_split(index_split, label_sets):
    """
    This is a helper function which creates a label split from an index split.
    @param index_split: a split of indices of the label sets
    @param label_sets: a set of labels for each index
    @return: a label split defined by the index split of the label sets
    """
    label_split = set()
    for index_selection in index_split:
        labels = set()
        for i in index_selection:
            labels.update(label_sets[i])
        label_split.add(frozenset(labels))
    return frozenset(label_split)

def split_svd(L_sqrt, epsilon=1e-14):
    """
    If a degenerate split is found then a DegenerateSplitException is raised.
    @param L_sqrt: the matrix square root of a Laplacian
    @param epsilon: small eigenvector loadings will be treated as zero
    @return: a set of two index sets defining a split of the indices
    """
    # get the fiedler vector
    v = sqrt_laplacian_to_fiedler(L_sqrt)
    # get the eigensplit
    eigensplit = eigenvector_to_split(v, epsilon)
    # validate the split
    min_cardinality, min_set = min((len(s), s) for s in eigensplit)
    if min_cardinality == 0:
        raise InvalidSpectralSplitException(D)
    elif min_cardinality == 1:
        index, = min_set
        raise DegenerateSplitException(index)
    else:
        return eigensplit

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


def sqrt_laplacian_to_fiedler(L_sqrt):
    """
    @param L_sqrt: the matrix square root of a Laplacian
    @return: the Fiedler vector of a related graph
    """
    U, S_array, VT = np.linalg.svd(L_sqrt, full_matrices=0)
    # we are interested in a column vector of U
    return U.T[-2]


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
        except DegenerateSplitException, e:
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
        next_ordered_labels = [ordered_labels[i] for i in index_selection]
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

def validate_data(data):
    """
    Each row of the input data is a gene.
    Each column of the input data is a genetic line of flies or something.
    @param data: row major list of lists of floating point numbers
    """
    if len(data) < 3:
        raise ValueError('too few rows of data')
    first_row_length = len(data[0])
    if first_row_length < 3:
        raise ValueError('too few columns of data')
    for row in data:
        if len(row) != first_row_length:
            raise ValueError('each row of data should have the same number of columns')

def get_data(filename):
    fin = open(filename)
    lines = fin.readlines()
    fin.close()
    arr = []
    for line in lines:
        line = line.strip()
        if not line:
            continue
        row = [float(x) for x in line.split()]
        arr.append(row)
    validate_data(arr)
    return np.array(arr)


"""
    def test_cool(self):
        X = get_data('cool.txt')
        L_sqrt = khorr.data_to_laplacian_sqrt(X)
        splits = get_splits(L_sqrt, split_using_svd, update_using_laplacian)
        #for s in splits:
            #print split_to_string(s)
"""

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
        X = get_data('kh-dataset.txt')
        L_sqrt = khorr.data_to_laplacian_sqrt(X)
        tree_data = TreeData(split_svd, update_svd)
        root = build_tree(L_sqrt, range(len(X)), tree_data)
        print
        print root.get_newick_string()
        print


if __name__ == '__main__':
    unittest.main()
