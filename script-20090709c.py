"""
This is a junk script.

Go back and analyze all of Eric's files that he did with MMC.
"""

import random
import os
import logging

import numpy as np

from khatrisvd import util
from khatrisvd import treebuilder
from khatrisvd import splitbuilder
from khatrisvd import khorr
from khatrisvd import heatmap
from khatrisvd import gradient

logging.basicConfig(level=logging.DEBUG)

g_output_directory = 'analysis-of-mmc-data-files'
g_input_directory = 'mmc-data-files'

def analyze(input_data_path, output_image_path, tree_building_function):
    logging.debug('read the matrix')
    X = util.file_to_comma_separated_matrix(input_data_path, has_headers=True)
    logging.debug('build the tree')
    root = tree_building_function(X)
    logging.debug('extract ordered indices from the tree')
    ordered_indices = root.ordered_labels()
    logging.debug('create the elementwise squared correlation matrix')
    R = np.corrcoef(X)
    logging.debug('permute the elementwise squared correlation matrix according to the ordering')
    M = heatmap.get_permuted_rows_and_columns(R, ordered_indices)
    logging.debug('create the heatmap')
    f = gradient.correlation_to_rgb
    heatmap.get_heatmap_with_dendrogram(M, root, f, output_image_path)

def main():
    for filename in os.listdir(g_input_directory):
        #if True:
        if filename.endswith('csv') and not filename.startswith('LGE'):
        #if filename.startswith('LGE'):
        #if filename.startswith('Sta'):
        #if filename.startswith('LGE') and filename.endswith('60.csv'):
            logging.debug(filename)
            input_data_path = os.path.join(g_input_directory, filename)
            # do the full analysis
            output_image_path = os.path.join(g_output_directory, filename + '.png')
            logging.debug('build the outgrouped tree')
            analyze(input_data_path, output_image_path, treebuilder.build_tree)
            # do a superficial analysis
            output_image_path = os.path.join(g_output_directory, filename + '.superficial.png')
            logging.debug('build a superficial tree from a single split')
            analyze(input_data_path, output_image_path, treebuilder.build_single_split_tree)
            # do a superficial analysis of correlation rather than squared correlation
            output_image_path = os.path.join(g_output_directory, filename + '.superficial.corr.png')
            logging.debug('build a superficial tree from a single split')
            analyze(input_data_path, output_image_path, treebuilder.build_single_split_correlation_tree)
            print

if __name__ == '__main__':
    main()
