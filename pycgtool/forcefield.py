import os

try:
    raise FileExistsError
except FileExistsError:
    pass
except NameError:
    class FileExistsError(Exception):
        pass


class ForceField:
    def __init__(self, name):
        try:
            os.makedirs(name)
        except FileExistsError as e:
            if not os.path.isdir(name):
                raise e

    def write_rtp(self, mapping, bonds):
        pass

