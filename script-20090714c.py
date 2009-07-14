"""
This is a junk script.

My code that directly constructs a downsampled heatmap is slow.
This script is written as an investigation of this phenomenon.
"""

import profile

import numpy as np

from khatrisvd import heatmap
from khatrisvd import gradient
from khatrisvd import khorr


def do_slow_stuff():
    # construct some random data
    p, n = 200, 10
    # construct the downsampled heatmap of the correlation of the data
    X = np.random.random((p, n))
    Z = khorr.get_standardized_matrix(X)
    color_function = gradient.correlation_to_rgb
    im = heatmap.get_reduced_heatmap_image(Z, color_function, reduction=10)

def main():
    profile.run('do_slow_stuff()')

if __name__ == '__main__':
    main()
