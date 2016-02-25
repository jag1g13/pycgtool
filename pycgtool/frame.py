"""
This module contains classes for storing atomic data.

The Frame class may contain multiple Residues which may each contain multiple Atoms.
Both Frame and Residue are iterable. Residue is indexable with either atom numbers or names.
"""

import numpy as np

from simpletraj import trajectory

np.seterr(all="raise")


class Atom:
    """
    Hold data for a single atom
    """
    __slots__ = ["name", "num", "type", "mass", "coords"]

    def __init__(self, name=None, num=None, type=None, mass=None, coords=None):
        self.name = name
        self.num = num
        self.type = type
        self.mass = mass
        self.coords = coords


class Residue:
    """
    Hold data for a residue - list of atoms
    """
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
        """
        Add an Atom to this Residue and store location in index

        :param atom: Atom to add to Residue
        :return: None
        """
        self.atoms.append(atom)
        self.name_to_num[atom.name] = len(self.atoms) - 1


class Frame:
    """
    Hold Atom data separated into Residues
    """

    def __init__(self, gro=None, xtc=None):
        """
        Return Frame instance having read Residues and Atoms from GRO if provided

        :param gro: GROMACS GRO file to read initial frame and extract residues
        :param xtc: GROMACS XTC file to read subsequent frames
        :return: Frame instance
        """
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
        """
        Read next frame from XTC

        :return: True if successful else False
        """
        try:
            self.xtc.get_frame(self.number)
            i = 0
            # Do this first outside the loop - numpy is fast
            self.xtc.x /= 10.
            x = self.xtc.x
            for res in self.residues:
                for atom in res:
                    atom.coords = x[i]
                    i += 1
            self.number += 1
            return True
        # IndexError - run out of xtc frames
        # AttributeError - we didn't provide an xtc
        except (IndexError, AttributeError):
            return False

    def _parse_gro(self, filename):
        """
        Parse a GROMACS GRO file and create Residues/Atoms
        Required before reading coordinates from XTC file

        :param filename: Filename of GROMACS GRO to read
        """
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
            self.number += 1

    def output_gro(self, filename):
        """
        Create a GROMACS GRO file from the data in this Frame

        :param filename: Name of GRO file to create
        """
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
        """
        Add a Residue to this Frame

        :param residue: Residue to add
        """
        self.residues.append(residue)
