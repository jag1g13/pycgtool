"""
This module contains classes for storing atomic data.

The Frame class may contain multiple Residues which may each contain multiple Atoms.
Both Frame and Residue are iterable. Residue is indexable with either atom numbers or names.
"""

import logging
import typing

import mdtraj
import numpy as np

from .parsers.cfg import CFG
from .util import backup_file, file_write_lines

logger = logging.getLogger(__name__)

np.seterr(all="raise")


class UnsupportedFormatException(Exception):
    def __init__(self, msg=None):
        if msg is None:
            msg = "Topology/Trajectory format not supported by this reader"
        super(UnsupportedFormatException, self).__init__(msg)


class NonMatchingSystemError(ValueError):
    def __init__(self, msg=None):
        if msg is None:
            msg = "Number of atoms does not match between topology and trajectory files"
        super(NonMatchingSystemError, self).__init__(msg)


class Trajectory:
    def __init__(self,
                 topology_file: typing.Optional[str] = None,
                 trajectory_file: typing.Optional[str] = None,
                 frame_start: int = 0):
        self.frame_number = frame_start

        if topology_file:
            try:
                self._trajectory = mdtraj.load(topology_file)
                self._topology = self._trajectory.topology

                if trajectory_file:
                    try:
                        self._trajectory = mdtraj.load(trajectory_file, top=self._topology)
                    
                    except ValueError as exc:
                        raise NonMatchingSystemError from exc

                self._load_coords(self.frame_number)
            
            except OSError as exc:
                if 'no loader' in str(exc):
                    raise UnsupportedFormatException from exc
                
                raise
        
        else:
            self._topology = mdtraj.Topology()

    def _load_coords(self, frame_number) -> None:
        for atom in self._topology.atoms:
            atom.coords = self._trajectory.xyz[frame_number, atom.index]

    def next_frame(self) -> bool:
        try:
            self._load_coords(self.frame_number + 1)
        
        except IndexError:
            return False

        self.frame_number += 1
        return True

    @property
    def residues(self):
        return self._topology.residues

    def residue(self, *args, **kwargs) -> mdtraj.core.topology.Residue:
        return self._topology.residue(*args, **kwargs)

    @property
    def atoms(self):
        return self._topology.atoms

    def atom(self, *args, **kwargs) -> mdtraj.core.topology.Atom:
        return self._topology.atom(*args, **kwargs)

    @property
    def natoms(self) -> int:
        return self._topology.n_atoms

    @property
    def numframes(self) -> int:
        # The MDTraj trajectory has the topology file as frame 0
        return self._trajectory.n_frames - 1

    def add_residue(self,
                    name,
                    chain: typing.Optional[mdtraj.core.topology.Residue] = None,
                    **kwargs) -> mdtraj.core.topology.Residue:
        if hasattr(self, '_trajectory'):
            raise TypeError('Cannot edit residues if a trajectory has been loaded')

        if chain is None:
            try:
                chain = self._topology.chain(0)
            
            except IndexError:
                chain = self._topology.add_chain()

        return self._topology.add_residue(name, chain, **kwargs)

    def add_atom(self, name: str,
                element: typing.Optional[mdtraj.element.Element],
                 residue: mdtraj.core.topology.Residue) -> mdtraj.core.topology.Atom:
        if hasattr(self, '_trajectory'):
            raise TypeError('Cannot edit atoms if a trajectory has been loaded')

        return self._topology.add_atom(name, element, residue)

    @property
    def box(self):
        return self._trajectory.unitcell_lengths[self.frame_number]

    def save(self, filename, **kwargs):
        self._trajectory.save(filename, **kwargs)

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
            from .framereader import get_frame_reader
            self._trajreader = get_frame_reader(gro, traj=xtc, frame_start=frame_start)

            self._trajreader.initialise_frame(self)

            if self._trajreader.num_atoms != self.natoms:
                raise AssertionError("Number of atoms does not match between gro and xtc files.")
            self.numframes += self._trajreader.num_frames

            if itp is not None:
                self._parse_itp(itp)


    @classmethod
    def instance_from_reader(cls, reader):
        """
        Return Frame instance initialised from existing FrameReader object

        :param FrameReader reader: FrameReader object
        :return: Frame instance
        """
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
                if "scipy" in repr(e):
                    e.msg = "XTC output with MDTraj also requires Scipy"
                else:
                    e.msg = "XTC output requires the module MDTraj (and probably Scipy)"
                raise

            backup_file(filename)
            self._xtc_buffer = mdtraj.formats.XTCTrajectoryFile(str(filename), mode="w")

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
        # TODO replace with actual ITP parser
        with CFG(filename) as itp:
            for molecule in itp:
                itpres = Residue(itp[molecule])
                for line in itp[molecule]["atoms"]:
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
        outputs = {"gro": self._get_gro_lines,
                   "lammps": self._get_lammps_data_lines}
        try:
            lines = outputs[format]()
            file_write_lines(filename, lines)
        except KeyError:
            print("ERROR: Invalid output format {0}, coordinates will not be output.".format(format))

    def _get_lammps_data_lines(self):
        """
        Return lines of LAMMPS DATA file.

        :return List[str]: Lines of DATA file containing current coordinates

        """
        raise NotImplementedError("LAMMPS Data output has not yet been implemented.")

    def _get_gro_lines(self):
        """
        Return lines of GRO file.

        :return List[str]: Lines of GRO file containing current coordinates
        """
        ret_lines = [
            self.name,
            "{0:5d}".format(self.natoms)
        ]

        i = 1
        format_string = "{0:5d}{1:5s}{2:>5s}{3:5d}{4:8.3f}{5:8.3f}{6:8.3f}"
        for res in self.residues:
            for atom in res:
                ret_lines.append(format_string.format(res.num, res.name, atom.name, i, *atom.coords))
                i += 1

        ret_lines.append("{0:10.5f}{1:10.5f}{2:10.5f}".format(*self.box))

        return ret_lines

    def add_residue(self, residue):
        """
        Add a Residue to this Frame

        :param residue: Residue to add
        """
        self.residues.append(residue)
