"""
This is a junk script.

The idea is to read some data and newick trees,
then give nested gene clusters within each data set,
then give clusters that are the same between data sets.
I hope that each experiment was done with the same set of genes.
The newick trees do not represent the phylogenetic relationship among the genes;
rather they represent a kind of inferred correlation structure.
"""

import sys
import logging

from khatrisvd import treebuilder
from khatrisvd import mtree

def get_gene_names(data_filename):
    """
    @param data_filename: the name of the data file
    @return: a list of gene names
    """
    # read the gene names from the file
    fin = open(data_filename)
    lines = fin.readlines()
    fin.close()
    # strip the whitespace from the ends of the lines
    lines = [line.strip() for line in lines]
    # remove empty lines
    lines = [line for line in lines if line]
    # remove the first line of headers
    lines = lines[1:]
    # split lines by commas
    rows = [[s.strip() for s in line.split(',')] for line in lines]
    # get only the first element of each line
    names = [row[0] for row in rows]
    # remove the quotes from the gene names
    names = [name[1:-1].strip() for name in names]
    # return the list
    return names

def summarize_data():
    """
    The first two sys.argv arguments are meant to be the two data files.
    """
    # get the first gene names
    data_filename = sys.argv[1]
    names_a = get_gene_names(data_filename)
    print len(names_a), 'names'
    # get the second gene names
    data_filename = sys.argv[2]
    names_b = get_gene_names(data_filename)
    print len(names_b), 'names'
    print
    # compare the names
    for a, b in zip(names_a, names_b)[:10]:
        print a
        print b
        print
    # compare all of the names
    nequal = 0
    ndifferent = 0
    for a, b in zip(names_a, names_b):
        if a == b:
            nequal += 1
        else:
            ndifferent += 1
    print nequal, 'equal'
    print ndifferent, 'different'

def get_tree(tree_filename):
    fin = open(tree_filename)
    multistring = fin.read()
    fin.close()
    return mtree.newick_to_mtree(multistring)

def analyze_clusters():
    """
    First arg is data file, second arg is newick file.
    """
    logging.debug('parse the gene names')
    data_filename = sys.argv[1]
    names = get_gene_names(data_filename)
    logging.debug('parse the tree')
    tree_filename = sys.argv[2]
    root = get_tree(tree_filename)
    logging.debug('center and sort the tree')
    root = treebuilder.center_and_sort_tree(root)
    logging.debug('get the clusters')
    clusters = treebuilder.get_clusters(root)
    logging.debug('show the corresponding genes')
    length_cluster_pairs = list(sorted((len(c), list(c)) for c in clusters))
    for l, c in length_cluster_pairs:
        #print len(c)
        print ', '.join(names[x] for x in c)

def analyze_consensus_clusters():
    logging.debug('parse the gene names')
    data_filename = sys.argv[1]
    names = get_gene_names(data_filename)
    logging.debug('get clusters from the first tree')
    tree_filename = sys.argv[2]
    root = get_tree(tree_filename)
    root = treebuilder.center_and_sort_tree(root)
    first_clusters = treebuilder.get_clusters(root)
    logging.debug('get clusters from the second tree')
    tree_filename = sys.argv[3]
    root = get_tree(tree_filename)
    root = treebuilder.center_and_sort_tree(root)
    second_clusters = treebuilder.get_clusters(root)
    logging.debug('get consensus clusters')
    consensus = set(frozenset(c) for c in first_clusters) & set(frozenset(c) for c in second_clusters)
    logging.debug('show the corresponding genes')
    length_cluster_pairs = list(sorted((len(c), list(c)) for c in consensus))
    for l, c in length_cluster_pairs:
        #print len(c)
        print ', '.join(names[x] for x in c)

def main():
    analyze_clusters()

if __name__ == '__main__':
    main()

