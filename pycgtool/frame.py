"""
This module contains classes for storing atomic data.

The Frame class may contain multiple Residues which may each contain multiple Atoms.
Both Frame and Residue are iterable. Residue is indexable with either atom numbers or names.
"""

import os
import abc
import logging
import itertools

import numpy as np

from .util import backup_file, FixedFormatUnpacker
from .parsers.cfg import CFG

logger = logging.getLogger(__name__)

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

    def __init__(self, name, num, type=None, mass=None, charge=None, coords=None):
        """
        Create an atom.
        
        :param str name: The name of the atom
        :param int num: The atom number
        :param str type: The atom type
        :param float mass: The mass of the atom
        :param float charge: The charge of the atom
        :param coords: The coordinates of the atom
        """
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

        for attr in ("type", "mass", "charge", "coords"):
            if getattr(self, attr) is None:
                setattr(self, attr, getattr(other, attr))


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


class FrameReader(metaclass=abc.ABCMeta):
    def __init__(self, topname, trajname=None, frame_start=0):
        self._topname = topname
        self._trajname = trajname
        self._frame_number = frame_start

        self.num_atoms = 0
        self.num_frames = 0

    def initialise_frame(self, frame):
        self._initialise_frame(frame)

    def read_next(self, frame):
        result = self.read_frame_number(self._frame_number, frame)
        if result:
            self._frame_number += 1
        return result

    def read_frame_number(self, number, frame):
        try:
            time, coords, box = self._read_frame_number(number)
            frame.time = time
            frame.box = box

            for atom, coord_line in zip(itertools.chain.from_iterable(frame.residues), coords):
                atom.coords = coord_line

        except (IndexError, AttributeError):
            # IndexError - run out of xtc frames
            # AttributeError - we didn't provide an xtc
            return False
        return True

    @abc.abstractmethod
    def _initialise_frame(self, frame):
        pass

    @abc.abstractmethod
    def _read_frame_number(self, number):
        pass


class FrameReaderSimpleTraj(FrameReader):
    def __init__(self, topname, trajname=None, frame_start=0):
        """
        Open input XTC file from which to read coordinates using simpletraj library.

        :param topname: MD topology file - not used
        :param trajname: MD trajectory file to read subsequent frames
        :param frame_start: Frame number to start on, default 0
        """
        FrameReader.__init__(self, topname, trajname, frame_start)

        from simpletraj import trajectory

        if trajname is not None:
            try:
                self._traj = trajectory.get_trajectory(trajname)
            except OSError as e:
                if not os.path.isfile(trajname):
                    raise FileNotFoundError(trajname) from e
                e.args = ("Error opening file '{0}'".format(trajname),)
                raise

            self.num_atoms = self._traj.numatoms
            self.num_frames = self._traj.numframes

    def _initialise_frame(self, frame):
        """
        Parse a GROMACS GRO file and create Residues/Atoms
        Required before reading coordinates from XTC file

        :param frame: Frame instance to initialise from GRO file
        """
        with open(self._topname) as gro:
            frame.name = gro.readline().strip()
            self.num_atoms = int(gro.readline())
            frame.natoms = self.num_atoms
            resnum_last = None
            atnum = 0

            unpacker = FixedFormatUnpacker("I5,2A5,5X,3F8",
                                           FixedFormatUnpacker.FormatStyle.Fortran)

            for _ in range(self.num_atoms):
                resnum, resname, atomname, *coords = unpacker.unpack(gro.readline())
                coords = np.array(coords, dtype=np.float32)

                if resnum != resnum_last:
                    frame.residues.append(Residue(name=resname,
                                                  num=resnum))
                    resnum_last = resnum
                    atnum = 0

                atom = Atom(name=atomname, num=atnum, coords=coords)
                frame.residues[-1].add_atom(atom)
                atnum += 1

            frame.box = np.array([float(x) for x in gro.readline().split()[0:3]], dtype=np.float32)

    def _read_frame_number(self, number):
        """
        Read next frame from XTC using simpletraj library.
        """
        self._traj.get_frame(number)
        # SimpleTraj uses Angstrom, we want nanometers
        xyz = self._traj.x / 10
        box = np.diag(self._traj.box)[0:3] / 10

        return self._traj.time, xyz, box


class FrameReaderMDTraj(FrameReader):
    def __init__(self, topname, trajname=None, frame_start=0):
        """
        Open input XTC file from which to read coordinates using mdtraj library.

        :param topname: GROMACS GRO file from which to read topology
        :param trajname: GROMACS XTC file to read subsequent frames
        """
        FrameReader.__init__(self, topname, trajname, frame_start)

        try:
            import mdtraj
        except ImportError as e:
            if "scipy" in e.msg:
                e.msg = "The MDTraj FrameReader also requires Scipy"
            else:
                e.msg = "The MDTraj FrameReader requires the module MDTraj (and probably Scipy)"
            raise
        logger.warning("WARNING: Using MDTraj which renames solvent molecules")

        try:
            if trajname is None:
                self._traj = mdtraj.load(topname)
            else:
                self._traj = mdtraj.load(trajname, top=topname)
        except OSError as e:
            if not os.path.isfile(topname):
                raise FileNotFoundError(topname) from e
            if not os.path.isfile(trajname):
                raise FileNotFoundError(trajname) from e
            e.args = ("Error opening file '{0}' or '{1}'".format(topname, trajname),)
            raise

        self.num_atoms = self._traj.n_atoms
        self.num_frames = self._traj.n_frames

    def _initialise_frame(self, frame):
        """
        Parse a GROMACS GRO file and create Residues/Atoms
        Required before reading coordinates from XTC file

        :param frame: Frame instance to initialise from GRO file
        """
        import mdtraj
        top = mdtraj.load(self._topname)

        frame.name = ""
        self.num_atoms = top.n_atoms
        frame.natoms = top.n_atoms

        frame.residues = [Residue(name=res.name, num=res.resSeq) for res in top.topology.residues]

        for atom in top.topology.atoms:
            new_atom = Atom(name=atom.name, num=atom.serial,
                            coords=top.xyz[0][atom.index])
            frame.residues[atom.residue.index].add_atom(new_atom)

        frame.box = top.unitcell_lengths[0]

    def _read_frame_number(self, number):
        """
        Read next frame from XTC using mdtraj library.
        """
        return self._traj.time[number], self._traj.xyz[number], self._traj.unitcell_lengths[number]


class FrameReaderMDAnalysis(FrameReader):
    def __init__(self, topname, trajname=None, frame_start=0):
        import MDAnalysis

        super().__init__(topname, trajname, frame_start)

        if trajname is None:
            self._traj = MDAnalysis.Universe(topname)
        else:
            self._traj = MDAnalysis.Universe(topname, trajname)

        self.num_atoms = self._traj.atoms.n_atoms
        self.num_frames = self._traj.trajectory.n_frames

    def _initialise_frame(self, frame):
        frame.name = ""
        frame.natoms = self.num_atoms

        import MDAnalysis
        topol = MDAnalysis.Universe(self._topname)
        frame.box = topol.dimensions[0:3] / 10.

        for res in topol.residues:
            residue = Residue(name=res.resname, num=res.resnum)
            for atom in res.atoms:
                residue.add_atom(Atom(name=atom.name, num=atom.id, coords=atom.position / 10.))
            frame.residues.append(residue)

    def _read_frame_number(self, number):
        traj_frame = self._traj.trajectory[number]
        return traj_frame.time, traj_frame.positions / 10., traj_frame.dimensions[0:3] / 10.


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
        self.name = ""
        self.residues = []
        self.number = frame_start - 1
        self.time = 0
        self.numframes = 0
        self.natoms = 0
        self.box = np.zeros(3, dtype=np.float32)

        self._xtc_buffer = None

        if gro is not None:
            open_xtc = {"simpletraj": FrameReaderSimpleTraj,
                        "mdtraj":     FrameReaderMDTraj}
            try:
                self._trajreader = open_xtc[xtc_reader](gro, xtc, frame_start=frame_start)
            except KeyError as e:
                e.args = ("XTC reader {0} is not a valid option.".format(xtc_reader))
                raise

            self._trajreader.initialise_frame(self)

            if self._trajreader.num_atoms != self.natoms:
                raise AssertionError("Number of atoms does not match between gro and xtc files.")
            self.numframes += self._trajreader.num_frames

            if itp is not None:
                self._parse_itp(itp)

    @classmethod
    def instance_from_reader(cls, reader):
        obj = cls()
        obj._trajreader = reader
        obj._trajreader.initialise_frame(obj)
        return obj

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
                atoms.append(repr(atom))
        rep += "\n".join(atoms)
        return rep

    def yield_resname_in(self, container):
        for res in self:
            if res.name in container:
                yield res

    def next_frame(self):
        """
        Read next frame from input XTC.

        :return: True if successful else False
        """
        result = self._trajreader.read_next(self)
        if result:
            self.number += 1
        return result

    def write_xtc(self, filename):
        """
        Write frame to output XTC file.

        :param filename: XTC filename to write to
        """
        if self._xtc_buffer is None:
            try:
                import mdtraj
            except ImportError as e:
                if "scipy" in e.msg:
                    e.msg = "XTC output with MDTraj also requires Scipy"
                else:
                    e.msg = "XTC output requires the module MDTraj (and probably Scipy)"
                raise

            backup_file(filename)
            self._xtc_buffer = mdtraj.formats.XTCTrajectoryFile(filename, mode="w")

        xyz = np.ndarray((1, self.natoms, 3), dtype=np.float32)
        i = 0
        for residue in self.residues:
            for atom in residue.atoms:
                xyz[0][i] = atom.coords
                i += 1

        time = np.array([self.time], dtype=np.float32)
        step = np.array([self.number], dtype=np.int32)

        box = np.zeros((1, 3, 3), dtype=np.float32)
        for i in range(3):
            box[0][i][i] = self.box[i]

        self._xtc_buffer.write(xyz, time=time, step=step, box=box)

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
        backup_file(filename)

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
