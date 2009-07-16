"""
This is supposed to be some GUI thing.

The Ubuntu package name for Image is python-imaging.
The Ubuntu package name for ImageTk is python-imaging-tk.
"""

import Tkinter

def main():
    root = Tkinter.Tk()
    container = Tkinter.Frame(root)
    container.pack()
    button = Tkinter.Button(container)
    button['text']= 'hello world'
    button['background'] = 'blue'
    button.pack()
    root.mainloop()

if __name__ == '__main__':
    main()
