"""
Make a progress bar.
"""

import os
import sys, time
from array import array
from fcntl import ioctl
import termios
import signal


class Bar:

    def __init__(self, high=100):
        """
        @param high: when the progress reaches this value then we are done
        """
        assert high > 0
        self.high = high
        self.outfile = sys.stderr
        self.finished = False
        # get the terminal width and reset the cached nfilled value
        self.on_resize(None, None)
        signal.signal(signal.SIGWINCH, self.on_resize)
        # set the progress to zero
        self.update(0)

    def on_resize(self, signal_number, frame):
        """
        @param signal_number: an unused parameter passed by the signal framework
        @param frame: an unused parameter passed by the signal framework
        """
        ioctl_result = ioctl(self.outfile, termios.TIOCGWINSZ, '\0'*8)
        nrows, self.ncols = array('h', ioctl_result)[:2]
        self.cached_nfilled = None

    def get_nfilled(self):
        """
        @return: the length of the filled bar
        """
        nbar_columns = self.ncols - 2
        nfilled = (self.progress * nbar_columns) / self.high
        return nfilled

    def get_progress_line(self):
        """
        @return: the progress string to print
        """
        nbar_columns = self.ncols - 2
        nfilled = self.get_nfilled()
        progress_line = '[' + '='*nfilled + '.'*(nbar_columns-nfilled) + ']'
        return progress_line

    def update(self, progress):
        """
        @param progress: the amount of progress made so far
        """
        assert 0 <= progress <= self.high
        if not self.finished:
            self.progress = progress
            nfilled = self.get_nfilled()
            if self.progress == self.high:
                self.cached_nfilled = nfilled
                progress_line = self.get_progress_line()
                self.outfile.write(progress_line + '\n')
                signal.signal(signal.SIGWINCH, signal.SIG_DFL)
                self.finished = True
            elif self.cached_nfilled != nfilled:
                self.cached_nfilled = nfilled
                progress_line = self.get_progress_line()
                self.outfile.write(progress_line + '\r')

    def finish(self):
        self.update(self.high)

