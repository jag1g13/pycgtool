import numpy as np

from simpletraj import trajectory

from .parsers.cfg import CFG


class BeadMap:
    __slots__ = ["name", "typ", "atoms"]

    def __init__(self, name=None, typ=None, atoms=None):
        self.name = name
        self.typ = typ
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
                    self._mappings[mol.name] = []
                    molmap = self._mappings[mol.name]
                    for name, typ, *atoms in mol:
                        newbead = BeadMap(name=name, typ=typ, atoms=atoms)
                        molmap.append(newbead)

    def __len__(self):
        return len(self._mappings)

    def __contains__(self, item):
        return item in self._mappings

    def __getitem__(self, item):
        return self._mappings[item]

    def __iter__(self):
        return iter(self._mappings)

    def apply(self, frame, exclude=set()):
        cgframe = Frame()
        cgframe.name = frame.name
        cgframe.box = frame.box
        cgframe.natoms = 0
        cgframe.residues = []

        for aares in frame:
            if aares.name in exclude:
                continue
            try:
                molmap = self._mappings[aares.name]
            except KeyError:
                print("Residue {0} not found in mapping - will be ignored.".format(aares.name))
                continue
            res = Residue(name=aares.name, num=aares.num)
            res.atoms = [Atom(name=atom.name, typ=atom.typ) for atom in molmap]

            # Perform mapping
            for i, (bead, map) in enumerate(zip(res, molmap)):
                res.name_to_num[bead.name] = i
                bead.coords = np.zeros(3)
                for atom in map:
                    bead.coords += aares[atom].coords
                bead.coords /= len(map)
                cgframe.natoms += 1

            cgframe.residues.append(res)
        return cgframe


class Bond:
    __slots__ = ["atoms", "values"]

    def __init__(self, atoms=None):
        self.atoms = atoms
        self.values = []


def angle(a, b, c=None):
    if c is None:
        c = np.cross(a, b)
    det = np.linalg.det([a, b, c])
    dot = np.dot(a, b)
    return np.arctan2(det, dot) - np.pi/2


class Measure:
    def __init__(self, filename):
        with CFG(filename) as cfg:
            self._molecules = {}
            for mol in cfg:
                self._molecules[mol.name] = []
                molbnds = self._molecules[mol.name]
                for atomlist in mol:
                    molbnds.append(Bond(atoms=atomlist))

    def apply(self, frame):
        def calc_length(res, atoms):
            vec = res[atoms[1]].coords - res[atoms[0]].coords
            return np.sqrt(np.sum(vec*vec))

        def calc_angle(res, atoms):
            # for atom in atoms:
                # print(res[atom].coords)
            veca = res[atoms[1]].coords - res[atoms[0]].coords
            vecb = res[atoms[2]].coords - res[atoms[1]].coords
            ret = np.degrees(angle(veca, vecb))
            # print(ret)
            return ret

        def calc_dihedral(res, atoms):
            raise NotImplementedError

        for res in frame:
            mol_meas = self._molecules[res.name]
            for bond in mol_meas:
                try:
                    calc = {2: calc_length,
                            3: calc_angle,
                            4: calc_dihedral}
                    val = calc[len(bond.atoms)](res, bond.atoms)
                    bond.values.append(val)
                except NotImplementedError:
                    pass

    def __len__(self):
        return len(self._molecules)

    def __contains__(self, item):
        return item in self._molecules

    def __getitem__(self, item):
        return self._molecules[item]

    def __iter__(self):
        return iter(self._molecules)


class Atom:
    __slots__ = ["name", "num", "typ", "coords"]

    def __init__(self, name=None, num=None, typ=None, coords=None):
        self.name = name
        self.num = num
        self.typ = typ
        self.coords = coords


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
    def __init__(self, gro=None, xtc=None):
        self.residues = []
        self.number = -1
        if gro is not None:
            self._parse_gro(gro)
        if xtc is not None:
            self.xtc = trajectory.XtcTrajectory(xtc)

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

    def next_frame(self):
        try:
            self.xtc.get_frame(self.number)
            i = 0
            x = self.xtc.x / 10.
            for res in self.residues:
                for atom in res:
                    atom.coords = x[i]
                    i += 1
            self.number += 1
            return True
        except IndexError:
            return False

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
