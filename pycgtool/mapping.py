import numpy as np

from .frame import Atom, Residue, Frame
from .parsers.cfg import CFG


class BeadMap:
    __slots__ = ["name", "typ", "type", "atoms", "charge"]

    def __init__(self, name=None, type=None, atoms=None, charge=0):
        self.name = name
        self.type = type
        self.atoms = atoms
        self.charge = charge

    def __iter__(self):
        return iter(self.atoms)

    def __len__(self):
        return len(self.atoms)


class EmptyBeadError(Exception):
    pass


class Mapping:
    def __init__(self, filename=None):
        self._mappings = {}

        if filename is not None:
            with CFG(filename) as cfg:
                for mol in cfg:
                    self._mappings[mol.name] = []
                    molmap = self._mappings[mol.name]
                    for name, typ, *atoms in mol:
                        newbead = BeadMap(name=name, type=typ, atoms=atoms)
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
            res.atoms = [Atom(name=atom.name, type=atom.type) for atom in molmap]

            # Perform mapping
            for i, (bead, bmap) in enumerate(zip(res, molmap)):
                res.name_to_num[bead.name] = i
                bead.coords = np.zeros(3)
                n = 0
                for atom in bmap:
                    try:
                        bead.coords += aares[atom].coords
                        n += 1
                    except KeyError:
                        # Atom does not exist in residue
                        pass
                # bead.coords = functools.reduce(lambda a,b: a + aares[b].coords, bmap, 0)
                try:
                    bead.coords /= n
                except FloatingPointError:
                    raise EmptyBeadError("Bead {0} in molecule {1} contains no atoms.".format(
                        bead.name, aares.name
                    ))
                cgframe.natoms += 1

            cgframe.residues.append(res)
        return cgframe

