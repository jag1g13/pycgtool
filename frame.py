
class Mapping:
    def __init__(self):
        pass


def mappings_from_file(file):
    mappings = {}
    with open(file) as f:
        pass


class Atom:
    def __init__(self, name=None, num=None, type=None, coords=None):
        self.name = name
        self.num = num
        self.type = type
        self.coords = coords


class Bead(Atom):
    def __init__(self, *args, **kwargs):
        super(Bead, self).__init__(**kwargs)


class Residue:
    def __init__(self, name=None, num=None):
        self.atoms = []
        self.name = name
        self.num = num

    def add_atom(self, atom):
        self.atoms.append(atom)


class Frame:
    def __init__(self, gro=None):
        self.residues = []
        if gro is not None:
            self._parse_gro(gro)

    def _parse_gro(self, gro):
        with open(gro) as g:
            self.name = g.readline()
            natoms = int(g.readline())
            i = 0
            resnum_last = -1
            for line in g:
                resnum = int(line[0:5])
                resname = line[5:10].strip()
                atomname = line[10:15].strip()
                coords = float(line[20:28]), float(line[28:36]), float(line[36:44])
                if resnum != resnum_last:
                    self.residues.append(Residue(name=resname,
                                                 num=resnum))
                    resnum_last = resnum
                    atnum = 0
                atom = Atom(name=atomname, num=atnum, coords=coords)
                self.residues[-1].add_atom(atom)
                if i >= natoms - 1:
                    break
                i += 1
                atnum += 1

    def add_residue(self, residue):
        self.residues.append(residue)