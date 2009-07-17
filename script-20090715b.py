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

Abbreviations are im for PIL image and tkim for a tk image object.
"""


import Tkinter
import sys

import ImageTk
import Image
import numpy as np

from khatrisvd import util
from khatrisvd import khorr
from khatrisvd import mtree
from khatrisvd import dendro


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

class GeneNames:
    """
    Show some gene names in a window.
    """

    def __init__(self, app, parent, ordered_names):
        """
        @param app: the application object
        @param parent: the parent container
        @param ordered_names: ordered gene names
        """
        self.app = app
        self.parent = parent
        self.ordered_names = ordered_names
        # make the canvas
        self.canvas = Tkinter.Canvas(self.parent, width=10, height=self.npixels)
        self.tkim = None
        #TODO stuff here

    def on_selection(self, index_range):
        begin, end = index_range
        selected_names = [self.ordered_names[i] for i in range(begin, end)]
        # make a container with a bunch of labels with the names
        self.container = Tkinter.Frame()
        labels = None
        #TODO stuff here
        self.app.repack()


class EastDendrogramWindow:
    """
    Draw a dendrogram to the right of a correlation window.
    """

    def __init__(self, app, parent, npixels, root, label_to_leaf):
        """
        @param app: the application object
        @param parent: the parent container
        @param npixels: the height of the canvas
        @param root: the root of the mtree
        @param label_to_leaf: a map from labels to leaves
        """
        self.app = app
        self.parent = parent
        self.npixels = npixels
        self.root = root
        self.label_to_leaf = label_to_leaf
        # make the canvas
        self.canvas = Tkinter.Canvas(self.parent, width=10, height=self.npixels)
        # initialize some junk that will be created when there is actually a dendrogram
        self.im = None
        self.nleaves = None
        self.tkim = None

    def on_selection(self, index_range):
        """
        Draw a new dendrogram.
        @param row_range: the begin and end indices
        """
        # create the stable subtree
        begin, end = index_range
        leaves = [self.label_to_leaf[label] for label in range(begin, end)]
        cloned = mtree.leaves_to_subtree(self.root, leaves)
        # create the dendrogram PIL image
        self.nleaves = len(leaves)
        blocksize = self.npixels / self.nleaves
        breadth_gap = blocksize - 1
        height_gap = 3
        image_height = self.npixels
        image_width = dendro.get_dendrogram_height(cloned, height_gap)
        print 'image width:', image_width
        print 'image height:', image_height
        self.im = Image.new('RGB', (image_width, image_height), 'white')
        dendro.draw_dendrogram(cloned, breadth_gap, height_gap, self.on_draw_line)
        # create the dendrogram tkinter image
        self.tkim = ImageTk.PhotoImage(self.im)
        # remake the canvas with the new image
        self.canvas.destroy()
        self.canvas = Tkinter.Canvas(self.parent, width=image_width, height=image_height)
        self.canvas.create_image(0, 0, image=self.tkim, anchor=Tkinter.NW)
        self.app.repack()

    def on_draw_line(self, line):
        """
        This is called by the dendrogram builder.
        This assumes that there is a self.im image of the correct size.
        @param line: each endpoint is a (breadth_offset, height_offset) pair
        """
        blocksize = self.npixels / self.nleaves
        initial_breadth_offset = blocksize / 2
        black = (0, 0, 0)
        ((a, b), (c, d)) = line
        if a == c:
            for height_offset in range(min(b, d), max(b, d) + 1):
                self.im.putpixel((height_offset, initial_breadth_offset + a), black)
        elif b == d:
            for breadth_offset in range(min(a, c), max(a, c) + 1):
                self.im.putpixel((b, initial_breadth_offset + breadth_offset), black)


class LowZoom:
    """
    Show the whole heatmap at a low zoom level.
    """

    def __init__(self, app, parent, npixels, image, nindices):
        """
        @param parent: the parent container
        @param npixels: the width and height of the canvas
        @param image: the low resolution image
        @param nindices: the number of rows and columns in the correlation matrix
        """
        # initialize the zoom target
        self.zoom_target_function = None
        # save the args
        self.app = app
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

    def __init__(self, app, parent, npixels, Z, nindices):
        """
        @param parent: the parent container
        @param npixels: the width and height of the canvas
        @param Z: None or this matrix times its transpose gives the correlation matrix
        @param nindices: the number of rows and column in the correlation matrix
        """
        self.app = app
        self.parent = parent
        self.npixels = npixels
        self.Z = Z
        self.nindices = nindices
        # do validation
        if self.Z:
            if len(self.Z) != self.nindices:
                raise ValueError('the number of rows in Z should be the same as the number of indices')
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
        low = index - self.npixels/2
        high = low + self.npixels
        if low < 0:
            low = 0
            high = low + self.npixels
        elif high > self.nindices:
            high = self.nindices
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
        if self.Z:
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

    def __init__(self, app, parent, npixels, Z, nindices):
        """
        @param parent: the parent container
        @param npixels: the width and height of the canvas
        @param Z: None or this matrix times its transpose gives the correlation matrix
        @param nindices: the number of rows and column in the correlation matrix
        """
        # save the args
        self.app = app
        self.parent = parent
        self.npixels = npixels
        self.Z = Z
        self.nindices = nindices
        # do validation
        if self.Z:
            if len(self.Z) != self.nindices:
                raise ValueError('the number of rows in Z should be the same as the number of indices')
        # each correlation coefficient is represented by a block this many pixels wide
        self.blocksize = 10
        # make the window
        self.canvas = Tkinter.Canvas(self.parent, width=self.npixels, height=self.npixels, bg='lightblue')
        # tell this function when a range of indices has been zoomed to
        self.selection_target_function = None

    def _center_to_range(self, index):
        """
        @param index: the center index
        @return: a (begin, end) range
        """
        # calculate the number of correlation coefficient rows and columns representable per screen
        nindices_visible = self.npixels / self.blocksize
        # get the range of visible indices
        low = index - nindices_visible / 2
        high = low + nindices_visible
        if low < 0:
            low = 0
            high = low + nindices_visible
        elif high > self.nindices:
            high = self.nindices
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
        if self.Z:
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
        # some indices have been selected
        if self.selection_target_function:
            self.selection_target_function(self.row_range, self.column_range)


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
        # create the list of ordered gene names
        print 'creating the list of ordered gene names...'
        self.mtree_root = mtree.newick_file_to_mtree(tree_filename)
        names = get_gene_names(data_filename)
        ordered_gene_names = [names[row_index] for row_index in self.mtree_root.ordered_labels()]
        # create the standardized and sorted data matrix
        #print 'creating the ordered correlation square root...'
        #X = util.file_to_comma_separated_matrix(data_filename, has_headers=True)
        #Z = khorr.get_standardized_matrix(X)
        #Z = np.vstack(Z[row_index] for row_index in self.mtree_root.ordered_labels())
        Z = None
        # create the label to leaf map
        print 'creating the map from labels to leaves...'
        self.label_to_leaf = dict((tip.label, tip) for tip in self.mtree_root.ordered_tips())
        # initialize the windows individually
        self._init_windows(low_zoom_image, Z, ordered_gene_names)
        # redo the layout of the windows
        self.repack()
        # initialize the connections among the windows
        self._connect_windows()
        # cache some of the information
        # TODO this information could go into other child windows
        self.ordered_gene_names = ordered_gene_names

    def _init_windows(self, low_zoom_image, Z, ordered_names):
        """
        Initialize the windows individually.
        @param low_zoom_image: a resized tkinter image
        @param Z: None or a standardized and sorted square root of the correlation matrix
        @param names: ordered gene names or probeset ids
        """
        nindices = len(ordered_names)
        self.low_zoom = LowZoom(self, self.parent, self.npixels, low_zoom_image, nindices)
        self.mid_zoom = MidZoom(self, self.parent, self.npixels, Z, nindices)
        self.high_zoom = HighZoom(self, self.parent, self.npixels, Z, nindices)
        self.east_dendrogram = EastDendrogramWindow(self, self.parent, self.npixels, self.mtree_root, self.label_to_leaf)

    def repack(self):
        """
        Redo the layouts.
        """
        self.low_zoom.canvas.pack(side=Tkinter.LEFT, fill=Tkinter.BOTH, expand=Tkinter.YES)
        self.mid_zoom.canvas.pack(side=Tkinter.LEFT, fill=Tkinter.BOTH, expand=Tkinter.YES)
        self.high_zoom.canvas.pack(side=Tkinter.LEFT, fill=Tkinter.BOTH, expand=Tkinter.YES)
        self.east_dendrogram.canvas.pack(side=Tkinter.LEFT, fill=Tkinter.BOTH, expand=Tkinter.YES)

    def _connect_windows(self):
        """
        Initialize the connections among the windows.
        """
        self.low_zoom.zoom_target_function = self.mid_zoom.on_zoom
        self.mid_zoom.zoom_target_function = self.high_zoom.on_zoom
        self.high_zoom.selection_target_function = self.on_selection

    def on_selection(self, row_index_range, column_index_range):
        """
        @param index_range: the begin and end indices of a range
        """
        self.east_dendrogram.on_selection(row_index_range)
        # FIXME this cheese_canvas stuff is temporary
        #self.cheese_canvas = Tkinter.Canvas(self.container, width=10, height=self.npixels, bg='lightblue')
        #self.cheese_canvas.pack(side=Tkinter.LEFT, fill=Tkinter.BOTH, expand=Tkinter.YES)
        # show stuff on the terminal
        print 'row genes:'
        begin, end = row_index_range
        for name in self.ordered_gene_names[begin:end]:
            print name
        print 'column genes:'
        begin, end = column_index_range
        for name in self.ordered_gene_names[begin:end]:
            print name


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