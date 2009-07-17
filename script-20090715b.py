"""
This will be a large matrix correlation viewer.

The Ubuntu package name for Image is python-imaging.
The Ubuntu package name for ImageTk is python-imaging-tk.

There should be three 300x300 pixel correlation windows laid out
horizontally, progressing from low zoom to high zoom.
The lowest zoom window should represent the entire matrix.
Each of the two higher zoom windows should each represent
some view of the matrix with full detail.
For the middle zoom window, each pixel is a correlation.
For the highest zoom window, each correlation gets a large
block of pixels roughly the size of a letter of text.
Gene names (or probeset ids or whatever) can be listed
next to rows of the highest zoom window.

Optionally the two high zoom windows could be accompanied by
dendrograms relating the corresponding columns of the correlation matrix.
A dendrogram for the low zoom window would be less useful.

The two low zoom windows could have crosshairs showing the region
represented at a higher zoom level in the window to its right.
Clicking on a point in a low zoom window should change the selected region.
"""


import Tkinter
import sys

import ImageTk
import Image
import numpy as np

from khatrisvd import util
from khatrisvd import khorr
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


class LowZoom:
    """
    Show the whole heatmap at a low zoom level.
    """

    def __init__(self, parent, npixels, image, nindices):
        """
        @param parent: the parent container
        @param npixels: the width and height of the canvas
        @param image: the low resolution image
        @param nindices: the number of rows and columns in the correlation matrix
        """
        # initialize the zoom target
        self.zoom_target_function = None
        # save the args
        self.parent = parent
        self.npixels = npixels
        self.image = image
        self.nindices = nindices
        # make the image
        self.canvas = Tkinter.Canvas(self.parent, width=self.npixels, height=self.npixels, bg='lightgreen')
        self.canvas.create_image(0, 0, image=image, anchor=Tkinter.NW)
        self.canvas.bind('<Button-1>', self.on_button_down)

    def on_button_down(self, event):
        """
        @param event: information about where the button was pressed
        """
        # get the point clicked on the window
        pt = (event.x, event.y)
        print 'clicked the low zoom window at', pt
        # if a target is registered then get the center indices
        if self.zoom_target_function:
            # calculate the center column from the x pixel
            center_column_index = (self.nindices * event.x) / self.npixels
            # calculate the center row from the y pixel
            center_row_index = (self.nindices * event.y) / self.npixels
            # call the zoom target
            self.zoom_target_function(center_row_index, center_column_index)


class MidZoom:
    """
    Show one pixel per correlation coefficient.
    """

    def __init__(self, parent, npixels, correlation_sqrt):
        """
        @param parent: the parent container
        @param npixels: the width and height of the canvas
        @param correlation_sqrt: this matrix times its transpose gives the correlation matrix
        """
        self.parent = parent
        self.npixels = npixels
        self.Z = correlation_sqrt
        self.canvas = Tkinter.Canvas(self.parent, width=self.npixels, height=self.npixels, bg='lightyellow')
        self.canvas.bind('<Button-1>', self.on_button_down)
        # initialize the image area
        self.mid_zoom_image = None
        self.row_range = None
        self.column_range = None
        # initialize the zoom target
        self.zoom_target_function = None

    def _center_to_range(self, index):
        """
        @param index: the center index
        @return: a (begin, end) range
        """
        nindices = len(self.Z)
        low = index - self.npixels/2
        high = low + self.npixels
        if low < 0:
            low = 0
            high = low + self.npixels
        elif high > nindices:
            high = nindices
            low = high - self.npixels
        return low, high

    def on_zoom(self, center_row_index, center_column_index):
        """
        @param center_row_index: the center row index if possible
        @param center_column_index: the center column index if possible
        """
        # get the new row and column ranges
        self.row_range = self._center_to_range(center_row_index)
        self.column_range = self._center_to_range(center_column_index)
        # convert the center indices to valid index ranges
        row_begin, row_end = self.row_range
        column_begin, column_end = self.column_range
        # create the correlation matrix
        R = np.dot(self.Z[row_begin:row_end], self.Z[column_begin:column_end].T)
        # create the red, green, and blue values
        red = np.minimum(np.maximum(1.5 - 2.0*np.abs(R - 0.5), 0.0), 1.0)*255
        green = np.minimum(np.maximum(1.5 - 2.0*np.abs(R), 0.0), 1.0)*255
        blue = np.minimum(np.maximum(1.5 - 2.0*np.abs(R + 0.5), 0.0), 1.0)*255
        # create the pil image
        im = Image.new('RGB', (self.npixels, self.npixels), 'green')
        for row_index in range(self.npixels):
            for column_index in range(self.npixels):
                offset = (row_index, column_index)
                rgb = (int(red[offset]), int(green[offset]), int(blue[offset]))
                im.putpixel((column_index, row_index), rgb)
        # create the tkinter image
        self.mid_zoom_image = ImageTk.PhotoImage(im)
        # draw the mid zoom image on the canvas
        self.canvas.create_image(0, 0, image=self.mid_zoom_image, anchor=Tkinter.NW)

    def on_button_down(self, event):
        """
        @param event: information about where the button was pressed
        """
        # get the point clicked on the window
        pt = (event.x, event.y)
        print 'clicked the mid zoom window at', pt
        # if a target is registered then get the center indices
        if self.zoom_target_function:
            if self.row_range and self.column_range:
                # unpack the row and column ranges
                column_begin, column_end = self.column_range
                row_begin, row_end = self.row_range
                # get the mouse click proportion
                column_proportion = event.x / float(self.npixels)
                row_proportion = event.y / float(self.npixels)
                # get the centers
                center_column_index = int(column_begin + self.npixels * column_proportion)
                center_row_index = int(row_begin + self.npixels * row_proportion)
                # call the zoom target
                self.zoom_target_function(center_row_index, center_column_index)


class HighZoom:
    """
    Show multiple pixels per correlation coefficient.
    """

    def __init__(self, parent, npixels, correlation_sqrt, ordered_names):
        """
        @param parent: the parent container
        @param npixels: the width and height of the canvas
        @param correlation_sqrt: this matrix times its transpose gives the correlation matrix
        @param ordered_names: ordered gene names
        """
        # save the args
        self.parent = parent
        self.npixels = npixels
        self.Z = correlation_sqrt
        self.ordered_names = ordered_names
        # each correlation coefficient is represented by a block this many pixels wide
        self.blocksize = 10
        # make the window
        self.canvas = Tkinter.Canvas(self.parent, width=self.npixels, height=self.npixels, bg='lightblue')

    def _center_to_range(self, index):
        """
        @param index: the center index
        @return: a (begin, end) range
        """
        # calculate the number of correlation coefficient rows and columns representable per screen
        nindices_visible = self.npixels / self.blocksize
        # calculate the total number of indices in the correlation matrix
        nindices_total = len(self.Z)
        # get the range of visible indices
        low = index - nindices_visible / 2
        high = low + nindices_visible
        if low < 0:
            low = 0
            high = low + nindices_visible
        elif high > nindices_total:
            high = nindices_total
            low = high - nindices_visible
        return low, high

    def on_zoom(self, center_row_index, center_column_index):
        """
        @param center_row_index: the center row index if possible
        @param center_column_index: the center column index if possible
        """
        # get the new row and column ranges
        self.row_range = self._center_to_range(center_row_index)
        self.column_range = self._center_to_range(center_column_index)
        # convert the center indices to valid index ranges
        row_begin, row_end = self.row_range
        column_begin, column_end = self.column_range
        # create the correlation matrix
        R = np.dot(self.Z[row_begin:row_end], self.Z[column_begin:column_end].T)
        # create the red, green, and blue values
        red = np.minimum(np.maximum(1.5 - 2.0*np.abs(R - 0.5), 0.0), 1.0)*255
        green = np.minimum(np.maximum(1.5 - 2.0*np.abs(R), 0.0), 1.0)*255
        blue = np.minimum(np.maximum(1.5 - 2.0*np.abs(R + 0.5), 0.0), 1.0)*255
        # calculate the number of correlation coefficient rows and columns representable per screen
        nindices_visible = self.npixels / self.blocksize
        # create the pil image
        im = Image.new('RGB', (self.npixels, self.npixels), 'green')
        for row_index in range(nindices_visible):
            for column_index in range(nindices_visible):
                offset = (row_index, column_index)
                rgb = (int(red[offset]), int(green[offset]), int(blue[offset]))
                for i in range(self.blocksize):
                    for j in range(self.blocksize):
                        x = column_index*self.blocksize + i
                        y = row_index*self.blocksize + j
                        im.putpixel((x, y), rgb)
        # create the tkinter image
        self.mid_zoom_image = ImageTk.PhotoImage(im)
        # draw the mid zoom image on the canvas
        self.canvas.create_image(0, 0, image=self.mid_zoom_image, anchor=Tkinter.NW)


class Main:

    def __init__(self, parent, data_filename, tree_filename, low_zoom_image_filename):
        """
        @param parent: the parent frame
        @param low_zoom_image: the path to a low resolution image
        """
        # save some args
        self.parent = parent
        # define some options
        self.npixels = 400
        # create the tkinter image
        print 'creating the low resolution image...'
        raw_pil_image = Image.open(low_zoom_image_filename)
        low_zoom_pil_image = raw_pil_image.resize((self.npixels, self.npixels), Image.ANTIALIAS)
        low_zoom_image = ImageTk.PhotoImage(low_zoom_pil_image)
        # create the standardized and sorted data matrix
        print 'creating the ordered correlation square root...'
        X = util.file_to_comma_separated_matrix(data_filename, has_headers=True)
        Z = khorr.get_standardized_matrix(X)
        root = mtree.newick_file_to_mtree(tree_filename)
        Z = np.vstack(Z[row_index] for row_index in root.ordered_labels())
        # create the list of ordered gene names
        print 'creating the list of ordered gene names...'
        names = get_gene_names(data_filename)
        ordered_gene_names = [names[row_index] for row_index in root.ordered_labels()]
        # make a container to help with the layout
        self.container = Tkinter.Frame(parent)
        self.container.pack()
        # initialize the windows individually
        self._init_windows(low_zoom_image, Z, ordered_gene_names)
        # initialize the connections among the windows
        self._connect_windows()

    def _init_windows(self, low_zoom_image, correlation_sqrt, ordered_names):
        """
        Initialize the windows individually.
        @param low_zoom_image: a resized tkinter image
        @param correlation_sqrt: a standardized and sorted square root of the correlation matrix
        @param names: ordered gene names or probeset ids
        """
        # initialize the low-zoom window
        nindices = len(ordered_names)
        self.low_zoom = LowZoom(self.container, self.npixels, low_zoom_image, nindices)
        self.low_zoom.canvas.pack(side=Tkinter.LEFT, fill=Tkinter.BOTH, expand=Tkinter.YES)
        # initialize the mid-zoom window
        self.mid_zoom = MidZoom(self.container, self.npixels, correlation_sqrt)
        self.mid_zoom.canvas.pack(side=Tkinter.LEFT, fill=Tkinter.BOTH, expand=Tkinter.YES)
        # initialize the high-zoom window
        self.high_zoom = HighZoom(self.container, self.npixels, correlation_sqrt, ordered_names)
        self.high_zoom.canvas.pack(side=Tkinter.LEFT, fill=Tkinter.BOTH, expand=Tkinter.YES)

    def _connect_windows(self):
        """
        Initialize the connections among the windows.
        """
        self.low_zoom.zoom_target_function = self.mid_zoom.on_zoom
        self.mid_zoom.zoom_target_function = self.high_zoom.on_zoom


def main():
    # get the filenames from the command line
    usage = 'python', sys.argv[0], '<data> <tree> <image>'
    if len(sys.argv) != 4:
        print usage
        return
    data_filename = sys.argv[1]
    tree_filename = sys.argv[2]
    low_zoom_image_filename = sys.argv[3]
    # initialize tkinter
    root = Tkinter.Tk()
    # create the gui and start the event loop
    main = Main(root, data_filename, tree_filename, low_zoom_image_filename)
    root.title('large correlation matrix viewer')
    print 'beginning the event loop'
    root.mainloop()

if __name__ == "__main__":
    main()
