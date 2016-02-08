

class Atom:
    def __init__(self, name=None, num=None, type=None):
        self.name = name
        self.num = num
        self.type = type

class Bead(Atom):
    def __init__(self, *args, **kwargs):
        super(Bead, self).__init__(**kwargs)

class Frame:
    def __init__(self):
        pass