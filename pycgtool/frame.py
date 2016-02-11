import numpy as np
from .parsers.cfg import CFG


class BeadMap:
    __slots__ = ["name", "type", "atoms"]

    def __init__(self, name=None, type=None, atoms=None):
        self.name = name
        self.type = type
        self.atoms = atoms

    def __iter__(self):
        return iter(self.atoms)

    def __len__(self):
        return len(self.atoms)


class Mapping:
    def __init__(self, filename=None):
        self._mappings = {}

        if filename is not None:
            with CFG(filename) as cfg:
                for mol in cfg:
                    if not mol.name in self._mappings:
                        self._mappings[mol.name] = []
                    for bead, type, *atoms in mol:
                        if bead not in self._mappings:
                            self._mappings[mol.name].append(BeadMap(name=bead, type=type, atoms=atoms))

    def __len__(self):
        return len(self._mappings)

    def __contains__(self, item):
        return item in self._mappings

    def __getitem__(self, item):
        return self._mappings[item]

    def apply(self, frame):
        cgframe = Frame()
        cgframe.name = frame.name
        cgframe.box = frame.box
        cgframe.natoms = 0
        cgframe.residues = [Residue(name=res.name, num=res.num) for res in frame.residues]
        for res, aares in zip(cgframe, frame):
            res.atoms = [Atom(name=atom.name, type=atom.type) for atom in self._mappings[res.name]]
            for bead, map in zip(res, self._mappings[res.name]):
                bead.coords = np.zeros(3)
                for atom in map:
                    bead.coords += aares[atom].coords
                bead.coords /= len(map)
                cgframe.natoms += 1
        return cgframe


class Atom:
    __slots__ = ["name", "num", "type", "coords"]

    def __init__(self, name=None, num=None, type=None, coords=None):
        self.name = name
        self.num = num
        self.type = type
        self.coords = coords


class Bead(Atom):
    def __init__(self, *args, **kwargs):
        super(Bead, self).__init__(**kwargs)


class Residue:
    __slots__ = ["name", "num", "atoms", "name_to_num"]

    def __init__(self, name=None, num=None):
        self.atoms = []
        self.name = name
        self.num = num
        self.name_to_num = {}

    def __iter__(self):
        return iter(self.atoms)

    def __getitem__(self, item):
        try:
            return self.atoms[item]
        except TypeError:
            return self.atoms[self.name_to_num[item]]

    def add_atom(self, atom):
        self.atoms.append(atom)
        self.name_to_num[atom.name] = len(self.atoms) - 1


class Frame:
    def __init__(self, gro=None):
        self.residues = []
        if gro is not None:
            self._parse_gro(gro)

    def __len__(self):
        return len(self.residues)

    def __iter__(self):
        return iter(self.residues)

    def __repr__(self):
        rep = self.name
        atoms = []
        for res in self.residues:
            for atom in res:
                atoms.append(repr(atom.coords))
        rep += "\n".join(atoms)
        return rep

    def _parse_gro(self, filename):
        with open(filename) as gro:
            self.name = gro.readline().strip()
            self.natoms = int(gro.readline())
            i = 0
            resnum_last = -1
            for line in gro:
                resnum = int(line[0:5])
                resname = line[5:10].strip()
                atomname = line[10:15].strip()
                coords = np.array([float(line[20:28]), float(line[28:36]), float(line[36:44])])
                if resnum != resnum_last:
                    self.residues.append(Residue(name=resname,
                                                 num=resnum))
                    resnum_last = resnum
                    atnum = 0
                atom = Atom(name=atomname, num=atnum, coords=coords)
                self.residues[-1].add_atom(atom)
                if i >= self.natoms - 1:
                    break
                i += 1
                atnum += 1
            line = gro.readline()
            self.box = np.array([float(x) for x in line.split()])

    def output_gro(self, filename):
        with open(filename, "w") as gro:
            print(self.name, file=gro)
            print("{0:5d}".format(self.natoms), file=gro)
            i = 1
            format_string = "{0:5d}{1:5s}{2:>5s}{3:5d}{4:8.3f}{5:8.3f}{6:8.3f}"
            for res in self.residues:
                for atom in res:
                    print(format_string.format(res.num, res.name,
                                               atom.name, i,
                                               *atom.coords), file=gro)
                    i += 1
            print("{0:10.5f}{1:10.5f}{2:10.5f}".format(*self.box), file=gro)

    def add_residue(self, residue):
        self.residues.append(residue)
