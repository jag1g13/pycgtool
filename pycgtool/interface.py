"""
This module contains classes for interaction at the terminal.
"""
import collections
import curses
import curses.textpad
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
        :param args: Optional program arguments from Argparse, will be displayed in interactive mode
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

    def interactive(self):
        """
        Read options in interactive terminal mode using curses.
        """
        curses.wrapper(self._inter)

    def _inter(self, stdscr):
        """
        Read options in interactive terminal mode using curses.

        :param stdscr: Curses window to use as interface
        """
        stdscr.clear()
        if self.args is not None:
            stdscr.addstr(1, 1, "Using GRO: {0}".format(self.args.gro))
            stdscr.addstr(2, 1, "Using XTC: {0}".format(self.args.xtc))
        stdscr.addstr(4, 1, "Press q to proceed")
        stdscr.box()
        stdscr.refresh()

        nrows = len(self)

        errscr = stdscr.derwin(3, curses.COLS - 3, nrows + 8, 1)
        errscr.border()

        window_config = stdscr.derwin(nrows + 2, curses.COLS - 3, 5, 1)
        window_config.box()
        window_config.refresh()
        window_keys = window_config.derwin(nrows, 20, 1, 0)
        window_config.vline(1, 18, curses.ACS_VLINE, nrows)
        window_vals = window_config.derwin(nrows, curses.COLS - 24, 1, 20)
        text_edit_wins = []
        text_inputs = []

        for i, (key, value) in enumerate(self):
            window_keys.addstr(i, 0, key)
            try:
                text_edit_wins.append(window_vals.derwin(1, 30, i, 0))
            except curses.error as e:
                raise RuntimeError("Your terminal is too small to fit the interface, please expand it") from e
            text_edit_wins[-1].addstr(0, 0, str(value))
            text_inputs.append(curses.textpad.Textbox(text_edit_wins[-1]))

        stdscr.refresh()
        window_keys.refresh()
        for window in text_edit_wins:
            window.refresh()

        pos = 0
        move = {"KEY_UP": lambda x: (x - 1) % nrows,
                "KEY_DOWN": lambda x: (x + 1) % nrows,
                "KEY_LEFT": lambda x: x,
                "KEY_RIGHT": lambda x: x}

        while True:
            key = text_edit_wins[pos].getkey(0, 0)
            errscr.erase()
            if key in move:
                pos = move[key](pos)
            if key == "\n":
                if type(self[pos]) is bool:
                    self._toggle_boolean_by_num(pos)
                else:
                    val = text_inputs[pos].edit().strip()
                    try:
                        self._set_by_num(pos, val)
                    except ValueError:
                        errscr.addstr(0, 0, "Invalid value '{0}' for option".format(val))
                        errscr.addstr(1, 0, "Value has been reset".format(val))

                text_edit_wins[pos].erase()
                text_edit_wins[pos].addstr(0, 0, str(self[pos]))
                text_edit_wins[pos].refresh()

            errscr.refresh()
            if key == "q":
                break


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
        self._start_time = time.clock()

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
        done = int(self._length * (self._its / self._maxits))
        left = self._length - done
        width = len(str(self._maxits))
        return "{0:-{width}} [".format(self._its, width=width) + done * "#" + left * "-" + "] {0}".format(self._maxits)

    def _stop(self):
        if not self._quiet:
            time_taken = int(time.clock() - self._start_time)
            print(self._bar + " took {0}s".format(time_taken))
        raise StopIteration

    def _display(self):
        try:
            time_remain = int((time.clock() - self._start_time) * ((self._maxits - self._its) / self._its))
        except ZeroDivisionError:
            time_remain = "-"
        print(self._bar + " {0}s left".format(time_remain), end="\r")
