"""This module contains classes for storing atomic data.

The Frame class may contain multiple Residues which may each contain multiple Atoms.
Both Frame and Residue are iterable. Residue is indexable with either atom numbers or names.
"""

import logging
import pathlib
import typing

import mdtraj
import numpy as np

logger = logging.getLogger(__name__)

np.seterr(all="raise")

PathLike = typing.Union[pathlib.Path, str]


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
                 topology_file: typing.Optional[PathLike] = None,
                 trajectory_file: typing.Optional[PathLike] = None,
                 frame_start: int = 0,
                 frame_end: typing.Optional[int] = None):
        """Load a simulation trajectory.

        :param topology_file: File containing system topology
        :param trajectory_file: File containing simulation trajectory
        :param frame_start: First frame of trajectory to use
        :param frame_end: Last frame of trajectory to use
        """
        if topology_file is not None:
            try:
                logging.info('Loading topology file')
                self._trajectory = mdtraj.load(str(topology_file))
                self._topology = self._trajectory.topology
                logging.info('Finished loading topology file')

                if trajectory_file is not None:
                    try:
                        logging.info('Loading trajectory file - this may take a while')
                        self._trajectory = mdtraj.load(str(trajectory_file),
                                                       top=self._topology)
                        self._trajectory = self._slice_trajectory(frame_start, frame_end)
                        logging.info(
                            'Finished loading trajectory file - loaded %d frames',
                            self._trajectory.n_frames)

                    except ValueError as exc:
                        raise NonMatchingSystemError from exc

                self._load_trajectory()

            except OSError as exc:
                if 'no loader' in str(exc):
                    raise UnsupportedFormatException from exc

                raise

        else:
            # No topology - we're probably building a CG frame
            self._topology = mdtraj.Topology()

    def _slice_trajectory(
            self,
            frame_start: int = 0,
            frame_end: typing.Optional[int] = None) -> mdtraj.Trajectory:
        """Slice input simulation trajectory to a subset of frames.

        :param frame_start: First frame of trajectory to use
        :param frame_end: Last frame of trajectory to use
        """
        traj = self._trajectory

        if frame_start != 0:
            if frame_end is not None:
                traj = traj[frame_start:frame_end]

            else:
                traj = traj[frame_start:]

        return traj

    def _load_trajectory(self) -> None:
        """Load a trajectory into the frame attributes."""
        # Improve performance by not double indexing repeatedly
        traj = self._trajectory

        for atom in self._topology.atoms:
            atom.coords = traj.xyz[:, atom.index]

        # TODO handle non-cubic boxes
        try:
            self.unitcell_lengths = traj.unitcell_lengths
            self.unitcell_angles = traj.unitcell_angles

        except TypeError:
            self.unitcell_lengths = None
            self.unitcell_angles = None

        self.time = traj.time

        # Cast unitcell vectors to float32 to avoid warning - they're currently float64
        if traj.unitcell_vectors is not None:
            traj.unitcell_vectors = traj.unitcell_vectors.astype(np.float32)

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
    def n_frames(self) -> int:
        """Number of frames in the trajectory."""
        return self._trajectory.n_frames

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

    def save(self,
             filename: PathLike,
             frame_number: typing.Optional[int] = None,
             **kwargs) -> None:
        """Write trajctory to file.

        :param filename: Name of output file
        """
        traj = self._trajectory
        if frame_number is not None:
            traj = traj.slice(frame_number)

        traj.save(str(filename), **kwargs)

    def build_trajectory(self) -> None:
        """Build an MDTraj trajectory from atom coordinates and the values stored as attributes on this frame."""
        xyz = np.array([atom.coords for atom in self._topology.atoms], dtype=np.float32)

        # We currently have axes: 0 - each atom, 1 - timestep, 2 - xyz coords
        # Need to convert to axes: 0 - timestep, 1 - each atom, 2 - xyz coords
        xyz = xyz.swapaxes(0, 1)

        optional_values = {
            attr: getattr(self, attr, None)
            for attr in {'time', 'unitcell_lengths', 'unitcell_angles'}
        }  # yapf: disable

        new_trajectory = mdtraj.Trajectory(xyz,
                                           topology=self._topology,
                                           **optional_values)

        if hasattr(self, '_trajectory'):
            self._trajectory += new_trajectory

        else:
            self._trajectory = new_trajectory
