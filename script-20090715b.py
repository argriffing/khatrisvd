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


import sys

import Tkinter
import ImageTk
import Image


class LowZoom:
    """
    Show the whole heatmap at a low zoom level.
    """

    def __init__(self, parent, npixels, image):
        """
        @param parent: the parent container
        @param npixels: the width and height of the canvas
        @param image: the low resolution image
        """
        # save the args
        self.parent = parent
        self.npixels = npixels
        self.image = image
        # make the image
        self.canvas = Tkinter.Canvas(self.parent, width=self.npixels, height=self.npixels, bg='lightgreen')
        self.canvas.create_image(0, 0, image=image, anchor=Tkinter.NW)
        self.canvas.bind('<Button-1>', self.on_button_down)

    def on_button_down(self, event):
        """
        @param event: information about where the button was pressed
        """
        print (event.x, event.y)


class MidZoom:
    """
    Show one pixel per correlation coefficient.
    """

    def __init__(self, parent, npixels):
        """
        @param parent: the parent container
        @param npixels: the width and height of the canvas
        """
        self.npixels = npixels
        self.parent = parent
        self.canvas = Tkinter.Canvas(self.parent, width=self.npixels, height=self.npixels, bg='lightyellow')
        self.canvas.bind('<Button-1>', self.on_button_down)

    def on_button_down(self, event):
        """
        @param event: information about where the button was pressed
        """
        print (event.x, event.y)


class HighZoom:
    """
    Show multiple pixels per correlation coefficient.
    """

    def __init__(self, parent, npixels):
        """
        @param parent: the parent container
        @param npixels: the width and height of the canvas
        """
        self.npixels = npixels
        self.parent = parent
        self.canvas = Tkinter.Canvas(self.parent, width=self.npixels, height=self.npixels, bg='lightblue')


class Main:

    def __init__(self, parent, low_zoom_image_filename):
        """
        @param parent: the parent frame
        @param low_zoom_image: the path to a low resolution image
        """
        # save some args
        self.parent = parent
        # define some options
        self.npixels = 400
        # create the tkinter image
        raw_pil_image = Image.open(low_zoom_image_filename)
        low_zoom_pil_image = raw_pil_image.resize((self.npixels, self.npixels), Image.ANTIALIAS)
        low_zoom_image = ImageTk.PhotoImage(low_zoom_pil_image)
        # make a container to help with the layout
        self.container = Tkinter.Frame(parent)
        self.container.pack()
        # put a canvas in the container
        self.low_zoom = LowZoom(self.container, self.npixels, low_zoom_image)
        self.low_zoom.canvas.pack(side=Tkinter.LEFT, fill=Tkinter.BOTH, expand=Tkinter.YES)
        # put a canvas in the container
        self.mid_zoom = MidZoom(self.container, self.npixels)
        self.mid_zoom.canvas.pack(side=Tkinter.LEFT, fill=Tkinter.BOTH, expand=Tkinter.YES)
        # put a canvas in the container
        self.high_zoom = HighZoom(self.container, self.npixels)
        self.high_zoom.canvas.pack(side=Tkinter.LEFT, fill=Tkinter.BOTH, expand=Tkinter.YES)


if __name__ == "__main__":
    # get the image filename from the command line
    low_zoom_image_filename = sys.argv[1]
    # initialize tkinter
    root = Tkinter.Tk()
    # create the gui and start the event loop
    main = Main(root, low_zoom_image_filename)
    root.title('large correlation matrix viewer')
    root.mainloop()

