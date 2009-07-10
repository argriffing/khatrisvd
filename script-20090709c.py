"""
This is a junk script.

Go back and analyze all of Eric's files that he did with MMC.
"""

import numpy as np
import random
import os

from khatrisvd import treebuilder
from khatrisvd import splitbuilder
from khatrisvd import khorr
from khatrisvd import heatmap
from khatrisvd import gradient

g_output_directory = 'analysis-of-mmc-data-files'
g_input_directory = 'mmc-data-files'

def data_file_to_matrix(input_data_path):
    """
    The input file is csv with headers.
    @param input_data_path: path to the input file
    @return: a matrix
    """
    # read the input file
    fin = open(input_data_path)
    lines = [line.strip() for line in fin.readlines()]
    fin.close()
    # skip empty lines
    lines = [line for line in lines if line]
    # skip the first line
    lines = lines[1:]
    # get rows of elements
    rows = [line.split(',') for line in lines]
    # skip the first element of each line
    rows = [row[1:] for row in rows]
    # convert elements to floats
    X = []
    for row_index, row in enumerate(rows):
        try:
            float_row = [float(x) for x in row]
        except ValueError, e:
            message_lines = [
                    'invalid number on row %d' % (row_index+1),
                    str(e)]
            raise ValueError('\n'.join(message_lines))
    rows = [[float(x) for x in row] for row in rows]
    # return the matrix
    return np.array(rows)

def analyze(input_data_path, output_image_path):
    print 'read the matrix'
    X = data_file_to_matrix(input_data_path)
    print 'get the initial sqrt of laplacian'
    L_sqrt = khorr.data_to_laplacian_sqrt(X)
    print 'create the tree'
    tree_data = treebuilder.TreeData(splitbuilder.split_svd, treebuilder.update_svd)
    root = treebuilder.build_tree(L_sqrt, range(len(X)), tree_data)
    print 'create the elementwise squared correlation matrix'
    R = np.corrcoef(X)
    #RoR = R*R
    print 'extract ordered indices from the tree'
    ordered_indices = root.ordered_labels()
    print 'permute the elementwise squared correlation matrix according to the ordering'
    #M = heatmap.get_permuted_rows_and_columns(RoR, ordered_indices)
    M = heatmap.get_permuted_rows_and_columns(R, ordered_indices)
    print 'create the heatmap'
    f = gradient.correlation_to_rgb
    heatmap.get_heatmap_with_dendrogram(M, root, f, output_image_path)

def main():
    for filename in os.listdir(g_input_directory):
        #if filename.endswith('csv') and not filename.startswith('LGE'):
        #if filename.startswith('LGE'):
        if filename.startswith('Sta'):
        #if filename.startswith('LGE') and filename.endswith('60.csv'):
            print filename
            input_data_path = os.path.join(g_input_directory, filename)
            output_image_path = os.path.join(g_output_directory, filename + '.png')
            analyze(input_data_path, output_image_path)
            print

if __name__ == '__main__':
    main()
