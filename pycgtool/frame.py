"""
This module contains classes for storing atomic data.

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
    def __init__(self, msg=None):
        if msg is None:
            msg = "Topology/Trajectory format not supported by this reader"
        super(UnsupportedFormatException, self).__init__(msg)


class NonMatchingSystemError(ValueError):
    def __init__(self, msg=None):
        if msg is None:
            msg = "Number of atoms does not match between topology and trajectory files"
        super(NonMatchingSystemError, self).__init__(msg)


class Frame:
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

                self._load_trajectory_frame(self.frame_number)

            except OSError as exc:
                if 'no loader' in str(exc):
                    raise UnsupportedFormatException from exc

                raise

        else:
            self._topology = mdtraj.Topology()

    def _load_trajectory_frame(self, frame_number) -> None:
        # Improve performance by not double indexing repeatedly
        xyz_frame = self._trajectory.xyz[frame_number]
        for atom in self._topology.atoms:
            atom.coords = xyz_frame[atom.index]

        # TODO handle non-cubic boxes
        try:
            self.unitcell_lengths = self._trajectory.unitcell_lengths[frame_number]
            self.unitcell_angles = self._trajectory.unitcell_angles[frame_number]

        except TypeError:
            self.unitcell_lengths = None
            self.unitcell_angles = None

        self.time = self._trajectory.time[frame_number]

    def next_frame(self) -> bool:
        try:
            self._load_trajectory_frame(self.frame_number + 1)

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

    def save(self, filename, **kwargs):
        self._trajectory.save(str(filename), **kwargs)

    def add_frame_to_trajectory(self) -> None:
        xyz = np.array([atom.coords for atom in self._topology.atoms])
        new_frame = mdtraj.Trajectory(xyz,
                                      topology=self._topology,
                                      time=[self.time],
                                      unitcell_lengths=[self.unitcell_lengths],
                                      unitcell_angles=[self.unitcell_angles])
        
        if hasattr(self, '_trajectory'):
            self._trajectory += new_frame
        
        else:
            self._trajectory = new_frame
