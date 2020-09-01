"""
This module contains classes for interaction at the terminal.
"""
import collections
import time


class Progress:
    """
    Display a progress bar during the main loop of a program.
    """

    def __init__(self, maxits, length=20, dowhile=None, quiet=False):
        """
        Return progress bar instance to handle printing of a progress bar within loops.

        :param maxits: Expected number of iterations
        :param length: Length of progress bar in characters
        :param dowhile: Function to call after each iteration, stops if False
        :param quiet: Skip printing of progress bar - for testing
        """
        self._maxits = maxits
        self._length = length
        self._dowhile = dowhile
        self._quiet = quiet
        self._its = -1
        self._start_time = time.process_time()

    def __len__(self):
        """
        Maximum iterator length.

        This length will be reached if the iterator is not stopped by a False dowhile condition or KeyboardInterrupt.

        :return: Maximum length of iterator
        """
        return self._maxits - self._its

    def __iter__(self):
        return self

    def __next__(self):
        """
        Allow iteration over Progress while testing dowhile condition.

        Will catch Ctrl-C and return control as if the iterator has been fully consumed.

        :return: Iteration number
        """
        self._its += 1

        try:
            if self._dowhile is not None and self._its > 0 and not self._dowhile():
                self._stop()
        except KeyboardInterrupt:
            print(end="\r")
            self._stop()

        if self._its >= self._maxits:
            self._stop()

        if not self._quiet:
            self._display()

        return self._its

    def run(self):
        """
        Iterate through self until stopped by maximum iterations or False condition.

        Use the tqdm library if it is present.
        """
        no_tqdm = False
        try:
            from tqdm import tqdm
        except ImportError:
            no_tqdm = True

        if self._quiet or no_tqdm:
            collections.deque(self, maxlen=0)

        else:
            self._quiet = True
            collections.deque(tqdm(self, total=len(self)-1, ncols=80), maxlen=0)
            self._quiet = False

        return self._its

    @property
    def _bar(self):
        done = self._length if self._maxits == 0 else int(self._length * (self._its / self._maxits))
        left = self._length - done
        width = len(str(self._maxits))
        return "{0:-{width}} [".format(self._its, width=width) + done * "#" + left * "-" + "] {0}".format(self._maxits)

    def _stop(self):
        if not self._quiet:
            time_elapsed = int(time.process_time() - self._start_time)
            print(self._bar + " took {0}s".format(time_elapsed))
            print(self._bar + " took {0} seconds".format(time_elapsed))

        raise StopIteration

    def _display(self):
        try:
            time_elapsed = (time.process_time() - self._start_time)
            proportion_remaining = (self._maxits - self._its) / self._its
            time_remain = int(time_elapsed * proportion_remaining)

        except ZeroDivisionError:
            time_remain = "-"
        print(self._bar + " {0}s left    ".format(time_remain), end="\r")
