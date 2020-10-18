"""This module contains classes for storing atomic data.

The Frame class may contain multiple Residues which may each contain multiple Atoms.
Both Frame and Residue are iterable. Residue is indexable with either atom numbers or names.
"""

import logging
import typing

import mdtraj
import numpy as np

logger = logging.getLogger(__name__)

np.seterr(all="raise")


class UnsupportedFormatException(Exception):
    """Exception raised when a topology/trajectory format cannot be parsed."""
    def __init__(self, msg=None):
        if msg is None:
            msg = "Topology/Trajectory format not supported by this reader"
        super(UnsupportedFormatException, self).__init__(msg)


class NonMatchingSystemError(ValueError):
    """Exception raised when topology and trajectory do not match."""
    def __init__(self, msg=None):
        if msg is None:
            msg = "Number of atoms does not match between topology and trajectory files"
        super(NonMatchingSystemError, self).__init__(msg)


class Frame:
    """Load and store data from a simulation trajectory."""
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
                        self._trajectory = mdtraj.load(trajectory_file,
                                                       top=self._topology)

                    except ValueError as exc:
                        raise NonMatchingSystemError from exc

                self._load_trajectory_frame(self.frame_number)

            except OSError as exc:
                if 'no loader' in str(exc):
                    raise UnsupportedFormatException from exc

                raise

        else:
            # No topology - we're probably building a CG frame
            self._topology = mdtraj.Topology()

    def _load_trajectory_frame(self, frame_number) -> None:
        """Load a trajectory frame into the frame attributes."""
        # Improve performance by not double indexing repeatedly
        traj = self._trajectory

        xyz_frame = traj.xyz[frame_number]
        for atom in self._topology.atoms:
            atom.coords = xyz_frame[atom.index]

        # TODO handle non-cubic boxes
        try:
            self.unitcell_lengths = traj.unitcell_lengths[frame_number]
            self.unitcell_angles = traj.unitcell_angles[frame_number]

        except TypeError:
            self.unitcell_lengths = None
            self.unitcell_angles = None

        self.time = traj.time[frame_number]

    def next_frame(self) -> bool:
        """Load the next trajectory frame into the frame attributes.

        :return: Loading next frame was successful?
        """
        try:
            self._load_trajectory_frame(self.frame_number + 1)

        except IndexError:
            return False

        self.frame_number += 1
        return True

    @property
    def residues(self):
        """Residues in the frame topology."""
        return self._topology.residues

    def residue(self, *args, **kwargs) -> mdtraj.core.topology.Residue:
        """Get a residue from the frame topology by name or id."""
        return self._topology.residue(*args, **kwargs)

    @property
    def atoms(self):
        """Atoms in the frame topology."""
        return self._topology.atoms

    def atom(self, *args, **kwargs) -> mdtraj.core.topology.Atom:
        """Get an atom from the frame topology by name or id."""
        return self._topology.atom(*args, **kwargs)

    @property
    def natoms(self) -> int:
        """Number of atoms in the frame topology."""
        return self._topology.n_atoms

    @property
    def numframes(self) -> int:
        """Number of frames in the frame trajectory."""
        # The MDTraj trajectory has the topology file as frame 0
        return self._trajectory.n_frames - 1

    def add_residue(self,
                    name,
                    chain: typing.Optional[
                        mdtraj.core.topology.Residue] = None,
                    **kwargs) -> mdtraj.core.topology.Residue:
        """Add a residue to the frame topology.

        :param name: Name of residue
        :param chain: MDTraj chain of residue
        :return: New residue
        """
        if hasattr(self, '_trajectory'):
            raise TypeError(
                'Cannot edit residues if a trajectory has been loaded')

        if chain is None:
            try:
                chain = self._topology.chain(0)

            except IndexError:
                chain = self._topology.add_chain()

        return self._topology.add_residue(name, chain, **kwargs)

    def add_atom(
            self, name: str, element: typing.Optional[mdtraj.element.Element],
            residue: mdtraj.core.topology.Residue
    ) -> mdtraj.core.topology.Atom:
        """Add an atom or CG bead to the frame topology.

        :param name: Name of atom
        :param element: MDTraj element of atom
        :param residue: MDTraj residue of atom
        :return: New atom
        """
        if hasattr(self, '_trajectory'):
            raise TypeError(
                'Cannot edit atoms if a trajectory has been loaded')

        return self._topology.add_atom(name, element, residue)

    def save(self, filename, **kwargs) -> None:
        """Write trajctory to file.

        :param filename: Name of output file
        """
        self._trajectory.save(str(filename), **kwargs)

    def add_frame_to_trajectory(self) -> None:
        """Add a new trajectory frame from the values stored as attributes on this frame."""
        xyz = np.array([atom.coords for atom in self._topology.atoms])

        optional_values = {
            attr: None if getattr(self, attr) is None else [getattr(self, attr)]
            for attr in {'time', 'unitcell_lengths', 'unitcell_angles'}
        }  # yapf: disable

        new_frame = mdtraj.Trajectory(xyz,
                                      topology=self._topology,
                                      **optional_values)

        if hasattr(self, '_trajectory'):
            self._trajectory += new_frame

        else:
            self._trajectory = new_frame
