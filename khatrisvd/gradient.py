"""
This has a gradient for the heatmap.
"""

import colorsys

g_correlation_gradient = [
        [0,0,0.5625],
        [0,0,0.625],
        [0,0,0.6875],
        [0,0,0.75],
        [0,0,0.8125],
        [0,0,0.875],
        [0,0,0.9375],
        [0,0,1],
        [0,0.0625,1],
		[0,0.125,1],
		[0,0.1875,1],
		[0,0.25,1],
		[0,0.3125,1],
		[0,0.375,1],
		[0,0.4375,1],
		[0,0.5,1],
		[0,0.5625,1],
		[0,0.625,1],
		[0,0.6875,1],
		[0,0.75,1],
		[0,0.8125,1],
		[0,0.875,1],
		[0,0.9375,1],
		[0,1,1],
		[0.0625,1,0.9375],
		[0.125,1,0.875],
		[0.1875,1,0.8125],
		[0.25,1,0.75],
		[0.3125,1,0.6875],
		[0.375,1,0.625],
		[0.4375,1,0.5625],
		[0.5,1,0.5],
		[0.5625,1,0.4375],
		[0.625,1,0.375],
		[0.6875,1,0.3125],
		[0.75,1,0.25],
		[0.8125,1,0.1875],
		[0.875,1,0.125],
		[0.9375,1,0.0625],
		[1,1,0],
		[1,0.9375,0],
		[1,0.875,0],
		[1,0.8125,0],
		[1,0.75,0],
		[1,0.6875,0],
		[1,0.625,0],
		[1,0.5625,0],
		[1,0.5,0],
		[1,0.4375,0],
		[1,0.375,0],
		[1,0.3125,0],
		[1,0.25,0],
		[1,0.1875,0],
		[1,0.125,0],
		[1,0.0625,0],
		[1,0,0],
		[0.9375,0,0],
		[0.875,0,0],
		[0.8125,0,0],
		[0.75,0,0],
		[0.6875,0,0],
		[0.625,0,0],
		[0.5625,0,0],
		[0.5,0,0]]

g_squared_correlation_gradient = []

def build_squared_correlation_gradient():
    global g_squared_correlation_gradient
    blue_hue, blue_saturation, value = colorsys.rgb_to_hsv(0, 0, 1.0)
    saturation = blue_saturation
    # create the image using M as the saturation value
    n = 1000
    for i in range(n):
        saturation = i / float(n - 1)
        colorsys_rgb = colorsys.hsv_to_rgb(blue_hue, saturation, value)
        pil_rgb = tuple(int(x*255) for x in colorsys_rgb)
        g_squared_correlation_gradient.append(pil_rgb)

def squared_correlation_to_rgb(rr):
    """
    @param r: a number between 0 and 1
    """
    nbins = len(g_squared_correlation_gradient)
    bin_index = int(rr * (nbins - 1) + 0.5)
    return tuple(g_squared_correlation_gradient[bin_index])

def correlation_to_rgb(r):
    """
    @param r: a number between -1 and 1
    """
    nbins = len(g_correlation_gradient)
    t = (r + 1.0) / 2.0
    bin_index = int(t * (nbins - 1) + 0.5)
    colorsys_rgb = g_correlation_gradient[bin_index]
    return tuple(int(x*255) for x in colorsys_rgb)

build_squared_correlation_gradient()
