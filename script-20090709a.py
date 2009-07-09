"""
This is a one-off script to check something Eric sent.
"""

from khatrisvd import mtree

def expand_internal_nodes(root):
    """
    Add a child node to each internal node with a label.
    The child node gets the label,
    and the label is removed from the internal node.
    """
    for node in root.preorder():
        if node.children and node.has_label():
            child = mtree.Node()
            child.label = node.label
            node.label = None
            node.add_child(child)

def doing_it_wrong_a():
    """
    In this function I tried to make a tree that I thought Eric meant.
    I made the wrong tree.
    """
    root = mtree.Node()
    last = root
    next_label = 0
    for i in range(32):
        for j in range(2):
            child = mtree.Node()
            child.label = next_label
            next_label += 1
            last.add_child(child)
        last = child
    expand_internal_nodes(root)
    print root.get_newick_string()

def doing_it_wrong_b():
    """
    In this function I tried to make a tree that I thought Eric meant.
    I made the wrong tree.
    """
    root = mtree.Node()
    shell = [root]
    next_label = 0
    for i in range(5):
        next_shell = []
        for node in shell:
            for j in range(2):
                child = mtree.Node()
                child.label = next_label
                next_label += 1
                node.add_child(child)
                next_shell.append(child)
        shell = next_shell
    expand_internal_nodes(root)
    print root.get_newick_string()

def main():
    root = mtree.Node()
    root.label = 0
    child = mtree.Node()
    child.label = 1
    root.add_child(child)
    nodes = [root, child]
    for i in range(1, 32):
        for j in (2*i, 2*i+1):
            child = mtree.Node()
            child.label = j
            nodes[i].add_child(child)
            nodes.append(child)
    expand_internal_nodes(root)
    print root.get_newick_string()

if __name__ == '__main__':
    main()
