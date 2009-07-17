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


class LowZoom:
    """
    Show the whole heatmap at a low zoom level.
    """

    def __init__(self, parent, npixels):
        """
        @param parent: the parent container
        @param npixels: the width and height of the canvas
        """
        self.npixels = npixels
        self.parent = parent
        self.canvas = Tkinter.Canvas(self.parent, width=self.npixels, height=self.npixels, bg='lightgreen')


class Main:

    def __init__(self, parent):
        # define some options
        self.npixels = 400
        # save the parent window
        self.parent = parent
        # make a container to help with the layout
        self.container = Tkinter.Frame(parent)
        self.container.pack()
        # put a canvas in the container
        self.low_zoom = LowZoom(self.container, self.npixels)
        self.low_zoom.canvas.pack(side=Tkinter.LEFT)
        # put a canvas in the container
        self.mid_zoom = LowZoom(self.container, self.npixels)
        self.mid_zoom.canvas.pack(side=Tkinter.LEFT)
        # put a canvas in the container
        self.high_zoom = LowZoom(self.container, self.npixels)
        self.high_zoom.canvas.pack(side=Tkinter.LEFT)

if __name__ == "__main__":
    root = Tkinter.Tk()
    main = Main(root)
    root.title('large correlation matrix viewer')
    root.mainloop()
