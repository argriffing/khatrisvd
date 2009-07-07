"""
This module is supposed to actually build the tree from the data.

Initially represent the tree as a set of informative splits.
Instead of repeatedly splitting a distance matrix,
the splitting function will repeatedly split a rectangular matrix.
"""

import mtree
import khorr

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

def split_using_svd(L_sqrt, epsilon=1e-14):
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

def update_using_laplacian(L_sqrt, index_set):
    """
    Update the matrix sqrt of Laplacian by summing rows of the removed indices.
    @param L_sqrt: the matrix square root of a Laplacian
    @param index_set: the set of indices that will be removed from the updated distance matrix
    @return: an updated matrix sqrt of Laplacian
    """
    return rmerge(L_sqrt, index_set)

def sqrt_laplacian_to_fiedler(L_sqrt):
    """
    @param L_sqrt: the matrix square root of a Laplacian
    @return: the Fiedler vector of a related graph
    """
    U, S_array, VT = np.linalg.svd(L_sqrt, full_matrices=0)
    # we are interested in a column vector of U
    return U.T[-2]

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

def get_splits(initial_L_sqrt, split_function, update_function, on_label_split=None):
    """
    Get the set of splits implied by the tree that would be reconstructed.
    @param initial_L_sqrt: the matrix square root of a Laplacian
    @param split_function: takes a matrix sqrt of Laplacian and returns an index split
    @param update_function: takes a matrix sqrt of Laplacian and an index subset and returns a matrix sqrt of Laplacian
    @param on_label_split: notifies the caller of the label split induced by an index split
    @return: a set of splits
    """
    n = len(initial_L_sqrt)
    # keep a stack of (label_set_per_vertex, distance_matrix) pairs
    initial_state = ([set([i]) for i in range(n)], initial_L_sqrt)
    stack = [initial_state]
    # process the stack in a depth first manner, building the split set
    label_split_set = set()
    while stack:
        label_sets, L_sqrt = stack.pop()
        # if the matrix is small then we are done
        if len(L_sqrt) < 4:
            continue
        # split the indices using the specified function
        try:
            index_split = split_function(L_sqrt)
            # convert the index split to a label split
            label_split = index_split_to_label_split(index_split, label_sets)
            # notify the caller if a callback is requested
            if on_label_split:
                on_label_split(label_split)
            # add the split to the master set of label splits
            label_split_set.add(label_split)
            # for large matrices create the new label sets and the new sqrt Laplacian matrices
            a, b = index_split
            for index_selection, index_complement in ((a,b),(b,a)):
                if len(index_complement) > 2:
                    next_label_sets = vmerge(label_sets, index_selection)
                    next_L_sqrt = update_function(L_sqrt, index_selection)
                    next_state = (next_label_sets, next_L_sqrt)
                    stack.append(next_state)
        except DegenerateSplitException, e:
            continue
    return label_split_set


class TestMe(unittest.TestCase):

    def assertAllClose(self, A, B, msg=None):
        """
        Assert that two matrices are almost entrywise identical.
        @param A: a matrix
        @param B: a matrix
        @param msg: shown when the assertion fails
        """
        self.assertTrue(np.allclose(A, B), msg=msg)

    def test_get_splits(self):
        p = 100
        n = 7
        X = np.random.random((p, n))
        L_sqrt = khorr.data_to_laplacian_sqrt(X)
        splits = get_splits(L_sqrt, split_using_svd, update_using_laplacian)
        for s in splits:
            print split_to_string(s)

    def test_rmerge(self):
        index_set = set([1, 2])
        M = np.array([[1,2,3,4],[5,6,7,8],[9,10,11,12],[13,14,15,16]])
        expected = np.array([[1,2,3,4],[14,16,18,20],[13,14,15,16]])
        observed = rmerge(M, index_set)
        self.assertAllClose(observed, expected)


if __name__ == '__main__':
    unittest.main()
