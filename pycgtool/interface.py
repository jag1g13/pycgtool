"""
This module contains classes for interaction at the terminal.
"""
import collections
import time


class Options:
    """
    Class to hold program options not specified at the initial command line.

    Values can be queried by indexing as a dictionary or by attribute.  Iterable.
    """
    def __init__(self, default, args=None):
        """
        Create Options instance from iterable of keys and default values.

        :param default: Iterable of key, default value pairs (e.g. list of tuples)
        :param args: Optional program arguments from Argparse
        """
        self._dict = collections.OrderedDict()
        for key, val in default:
            try:
                val = val.lower()
            except AttributeError:
                pass

            self._dict[key.lower()] = (val, type(val))

        # Allow to carry options from argparse
        self.args = args

    def __getattr__(self, attr):
        return self._dict[attr.lower()][0]

    def __repr__(self):
        res = "[" + ", ".join((str((key, val[0])) for key, val in self._dict.items())) + "]"
        return res

    def __iter__(self):
        return iter(((key, val[0]) for key, val in self._dict.items()))

    def __len__(self):
        return len(self._dict)

    def __getitem__(self, item):
        try:
            return self._dict[item]
        except KeyError:
            try:
                opt = list(self._dict.keys())[item]
                return self._dict[opt][0]
            except TypeError:
                raise TypeError("Must access Options using either a string or an integer")

    def set(self, opt, val):
        """
        Set an argument by name.

        :param opt: Option to set
        :param val: Value to set option to
        """
        opt = opt.lower()
        try:
            val = val.lower()
        except AttributeError:
            pass
        _type = self._dict[opt][1]

        if _type is not type(val):
            if _type is bool:
                self._dict[opt] = (_truthy(val), bool)
            else:
                self._dict[opt] = (_type(val), _type)
        else:
            self._dict[opt] = (val, _type)

    def _set_by_num(self, opt_num, val):
        """
        Set an argument if only its position in sequence is known.
        For use in Options._inter.

        :param opt_num: Sequence number of option to set
        :param val: Value to set option to
        """
        opt = list(self._dict.keys())[opt_num]
        self.set(opt, val)

    def toggle_boolean(self, opt):
        """
        Toggle a boolean argument by name.

        :param opt: Option to toggle
        """
        entry = self._dict[opt]
        if entry[1] is bool:
            self._dict[opt] = (not entry[0], entry[1])
        else:
            raise TypeError("Only boolean options can be toggled")

    def _toggle_boolean_by_num(self, opt_num):
        """
        Toggle a boolean argument if only its position in sequence is known.
        For use in Options._inter.

        :param opt_num: Sequence number of option to toggle
        """
        opt = list(self._dict.keys())[opt_num]
        self.toggle_boolean(opt)


def _truthy(string):
    """
    Evaluate a string as True or False in the natural way.

    :param string: String to evaluate
    :return: True or False
    """
    truthy_strings = ("yes", "y", "on", "true", "t", "1")
    falsey_strings = ("no", "n", "off", "false", "f", "0")

    string = string.lower().strip()
    if string in truthy_strings:
        return True
    elif string in falsey_strings:
        return False
    else:
        raise ValueError("Value '{0}' could not be converted to boolean".format(string))


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
