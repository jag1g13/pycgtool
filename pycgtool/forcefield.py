import os

# Python 3.2 doesn't have FileExistsError
try:
    raise FileExistsError
except NameError:
    class FileExistsError(Exception):
        pass
except FileExistsError:
    pass


class ForceField:
    def __init__(self, name):
        try:
            os.makedirs(name)
        except FileExistsError as e:
            if not os.path.isdir(name):
                raise e
        except OSError as e:
            if e.errno != 17 or not os.path.isdir(name):
                raise e

    def write_rtp(self, mapping, bonds):
        pass
