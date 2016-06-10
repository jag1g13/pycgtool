"""
This module contains classes for storing atomic data.

The Frame class may contain multiple Residues which may each contain multiple Atoms.
Both Frame and Residue are iterable. Residue is indexable with either atom numbers or names.
"""

import os

import numpy as np

from simpletraj import trajectory

from mdtraj.formats import XTCTrajectoryFile

from .util import backup_file
from .parsers.cfg import CFG

np.seterr(all="raise")


# Create FileNotFoundError if using older version of Python
try:
    try:
        raise FileNotFoundError
    except FileNotFoundError:
        pass
except NameError:
    class FileNotFoundError(OSError):
        pass


class Atom:
    """
    Hold data for a single atom
    """
    __slots__ = ["name", "num", "type", "mass", "charge", "coords"]

    def __init__(self, name=None, num=None, type=None, mass=None, charge=None, coords=None):
        self.name = name
        self.num = num
        self.type = type
        self.mass = mass
        self.charge = charge
        self.coords = coords

    def __repr__(self):
        return "Atom #{0} {1} type: {2} mass: {3} charge: {4}".format(
            self.num, self.name, self.type, self.mass, self.charge
        )

    def add_missing_data(self, other):
        assert self.name == other.name
        assert self.num == other.num
        if self.type is None:
            self.type = other.type
        if self.mass is None:
            self.mass = other.mass
        if self.charge is None:
            self.charge = other.charge
        if self.coords is None:
            self.coords = other.coords


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
            return self.atoms[self.name_to_num[item]]
        except KeyError:
            pass

        try:
            return self.atoms[item]
        except TypeError as e:
            e.args = ("Atom {0} does not exist in residue {1}".format(item, self.name),)
            raise

    def __len__(self):
        return len(self.atoms)

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

    def __init__(self, gro=None, xtc=None, itp=None, frame_start=0, xtc_reader="simpletraj"):
        """
        Return Frame instance having read Residues and Atoms from GRO if provided

        :param gro: GROMACS GRO file to read initial frame and extract residues
        :param xtc: GROMACS XTC file to read subsequent frames
        :param itp: GROMACS ITP file to read masses and charges
        :return: Frame instance
        """
        self.residues = []
        self.number = frame_start - 1
        self.numframes = 0
        self.box = np.zeros(3, dtype=np.float32)

        self._xtc_buffer = None

        if gro is not None:
            self._parse_gro(gro)
            self.numframes += 1

            if xtc is not None:
                self._xtc_reader = xtc_reader
                open_xtc = {"simpletraj": self._open_xtc_simpletraj,
                            "mdtraj":     self._open_xtc_mdtraj}
                try:
                    open_xtc[self._xtc_reader](xtc)
                except KeyError as e:
                    e.args = ("XTC reader {0} is not a valid option.".format(self._xtc_reader))
                    raise

            if itp is not None:
                self._parse_itp(itp)

    def _open_xtc_simpletraj(self, xtc):
        try:
            self.xtc = trajectory.XtcTrajectory(xtc)
        except OSError as e:
            if not os.path.isfile(xtc):
                raise FileNotFoundError(xtc) from e
            e.args = ("Error opening file '{0}'".format(xtc),)
            raise
        else:
            if self.xtc.numatoms != self.natoms:
                raise AssertionError("Number of atoms does not match between gro and xtc files.")
            self.numframes += self.xtc.numframes

    def _open_xtc_mdtraj(self, xtc):
        try:
            self.xtc = XTCTrajectoryFile(xtc)
        except OSError as e:
            if not os.path.isfile(xtc):
                raise FileNotFoundError(xtc) from e
            e.args = ("Error opening file '{0}'".format(xtc),)
            raise
        else:
            xyz, time, step, box = self.xtc.read(n_frames=1)
            natoms = len(xyz[0])
            if natoms != self.natoms:
                print(xyz[0])
                print(natoms, self.natoms)
                raise AssertionError("Number of atoms does not match between gro and xtc files.")

            # Seek to end to count frames
            # self.xtc.seek(0, whence=2)
            self.numframes += 1
            self.xtc.seek(0)
            while True:
                try:
                    self.xtc.seek(1, whence=1)
                    self.numframes += 1
                except IndexError:
                    break
            self.xtc.seek(0)

    def __len__(self):
        return len(self.residues)

    def __iter__(self):
        return iter(self.residues)

    def __getitem__(self, item):
        return self.residues[item]

    def __repr__(self):
        rep = self.name + "\n"
        atoms = []
        for res in self.residues:
            for atom in res:
                atoms.append(repr(atom.coords))
        rep += "\n".join(atoms)
        return rep

    def next_frame(self, exclude=None):
        """
        Read next frame from XTC

        :return: True if successful else False
        """
        xtc_read = {"simpletraj": self._next_frame_simpletraj,
                    "mdtraj":     self._next_frame_mdtraj}
        try:
            return xtc_read[self._xtc_reader](exclude)
        except KeyError as e:
            e.args = ("XTC reader {0} is not a valid option.".format(self._xtc_reader))
            raise

    def _next_frame_mdtraj(self, exclude=None):
        try:
            # self.xtc.seek(self.number)
            i = 0
            xyz, time, step, box = self.xtc.read(n_frames=1)
            xyz = xyz[0]
            for res in self.residues:
                if exclude is not None and res.name in exclude:
                    continue
                for atom in res:
                    atom.coords = xyz[i]
                    i += 1

            self.number += 1
            self.box = np.diag(box[0]) / 10.
            return True
        # IndexError - run out of xtc frames
        # AttributeError - we didn't provide an xtc
        except (IndexError, AttributeError):
            return False

    def _next_frame_simpletraj(self, exclude=None):
        try:
            self.xtc.get_frame(self.number)
            i = 0
            # Do this first outside the loop - numpy is fast
            self.xtc.x /= 10.
            for res in self.residues:
                if exclude is not None and res.name in exclude:
                    continue
                for atom in res:
                    atom.coords = self.xtc.x[i]
                    i += 1
            self.number += 1

            self.box = np.diag(self.xtc.box)[0:3] / 10

            return True
        # IndexError - run out of xtc frames
        # AttributeError - we didn't provide an xtc
        except (IndexError, AttributeError):
            return False

    class XTCBuffer:
        def __init__(self):
            self.coords = []
            self.step = []
            self.box = []

        def __call__(self):
            return (np.array(self.coords, dtype=np.float32),
                    np.array(self.step,   dtype=np.int32),
                    np.array(self.box,    dtype=np.float32))

        def append(self, coords, step, box):
            self.coords.append(coords)
            self.step.append(step)
            self.box.append(box)

    def write_to_xtc_buffer(self):
        if self._xtc_buffer is None:
            self._xtc_buffer = self.XTCBuffer()

        xyz = np.ndarray((self.natoms, 3))
        i = 0
        for residue in self.residues:
            for atom in residue.atoms:
                xyz[i] = atom.coords
                i += 1

        box = np.zeros((3, 3), dtype=np.float32)
        for i in range(3):
            box[i][i] = self.box[i]

        self._xtc_buffer.append(xyz, self.number, box)

    def flush_xtc_buffer(self, filename):
        if self._xtc_buffer is not None:
            xtc = XTCTrajectoryFile(filename, mode="w")

            xyz, step, box = self._xtc_buffer()
            xtc.write(xyz, step=step, box=box)
            xtc.close()

            self._xtc_buffer = None

    def write_xtc(self, filename):
        if self._xtc_buffer is None:
            self._xtc_buffer = XTCTrajectoryFile(filename, mode="w")

        xyz = np.ndarray((1, self.natoms, 3), dtype=np.float32)
        i = 0
        for residue in self.residues:
            for atom in residue.atoms:
                xyz[0][i] = atom.coords
                i += 1

        step = np.array([self.number], dtype=np.int32)

        box = np.zeros((1, 3, 3), dtype=np.float32)
        for i in range(3):
            box[0][i][i] = self.box[i]

        self._xtc_buffer.write(xyz, step=step, box=box)
        # self._xtc_buffer.close()

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
                coords = np.array([float(line[20:28]), float(line[28:36]), float(line[36:44])],
                                  dtype=np.float32)

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
            self.box = np.array([float(x) for x in line.split()[0:3]], dtype=np.float32)
            self.number += 1

    def _parse_itp(self, filename):
        """
        Parse a GROMACS ITP file to extract atom charges/masses.

        Optional but requires that ITP contains only a single residue.

        :param filename: Filename of GROMACS ITP to read
        """
        with CFG(filename) as itp:
            itpres = Residue(itp["moleculetype"][0][0])
            for line in itp["atoms"]:
                atom = Atom(num=int(line[0]) - 1, type=line[1], name=line[4], charge=float(line[6]), mass=float(line[7]))
                itpres.add_atom(atom)

            for res in self.residues:
                if res.name == itpres.name:
                    for atom, itpatom in zip(res, itpres):
                        atom.add_missing_data(itpatom)

    def output(self, filename, format="gro"):
        """
        Write coordinates from Frame to file.

        :param filename: Name of file to write to
        :param format: Format to write e.g. 'gro', 'lammps'
        """
        outputs = {"gro": self._output_gro,
                   "lammps": self._output_lammps_data}
        try:
            outputs[format](filename)
        except KeyError:
            print("ERROR: Invalid output format {0}, coordinates will not be output.".format(format))

    def _output_lammps_data(self, filename):
        """
        Output Frame coordinates in LAMMPS data format.

        :param filename: Name of DATA file to create
        """
        raise NotImplementedError("LAMMPS Data output has not yet been implemented.")

    def _output_gro(self, filename):
        """
        Create a GROMACS GRO file from the data in this Frame

        :param filename: Name of GRO file to create
        """
        backup_file(filename, verbose=True)

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
