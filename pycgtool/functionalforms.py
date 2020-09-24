"""Utilities for describing the functional forms of bonded potentials being calculated.

See http://manual.gromacs.org/documentation/current/reference-manual/functions/bonded-interactions.html
and http://manual.gromacs.org/documentation/current/reference-manual/topologies/topology-file-formats.html#tab-topfile2

The links above list the bonded potentials available for use in GROMACS force fields.
These can be implemented as functional forms in PyCGTOOL by subclassing :class:`FunctionalForm`.
"""
import abc
import enum
import math
import typing

import numpy as np


def get_functional_forms() -> enum.Enum:
    """Get enum of known functional forms."""
    enum_dict = {}
    for subclass in FunctionalForm.__subclasses__():
        name = subclass.__name__
        enum_dict[name] = subclass

    return enum.Enum('FunctionalForms', enum_dict)


class FunctionalForm(metaclass=abc.ABCMeta):
    """Parent class of any functional form used in Boltzmann Inversion to convert variance to a force constant."""
    def __init__(self,
                 mean_func: typing.Callable = np.nanmean,
                 variance_func: typing.Callable = np.nanvar):
        """Inject functions for calculating the mean and variance into the Boltzmann Inversion equations.

        :param mean_func: Function to calculate the mean - default is np.nanmean
        :param variance_func: Function to calculate the variance - default is np.nanvar
        """
        self._mean_func = mean_func
        self._variance_func = variance_func

    def eqm(self, values: typing.Iterable[float], temp: float) -> float:  # pylint: disable=unused-argument
        """Calculate equilibrium value.

        May be overridden by functional forms.

        :param values: Measured internal coordinate values from which to calculate equilibrium value
        :param temp: Temperature of simulation
        :return: Calculated equilibrium value
        """
        return self._mean_func(values)

    @abc.abstractmethod
    def fconst(self, values: typing.Iterable[float], temp: float) -> float:
        """Calculate force constant.
        Abstract static method to be defined by all functional forms.

        :param values: Measured internal coordinate values from which to calculate force constant
        :param temp: Temperature of simulation
        :return: Calculated force constant
        """
        raise NotImplementedError

    @abc.abstractproperty
    def gromacs_type_ids(self) -> typing.Tuple[int]:
        """Return tuple of GROMACS potential type ids when used as length, angle, dihedral.

        :return tuple[int]: Tuple of GROMACS potential type ids
        """
        raise NotImplementedError

    @classmethod
    def gromacs_type_id_by_natoms(cls, natoms: int) -> int:
        """Return the GROMACS potential type id for this functional form when used with natoms.

        :param int natoms: Number of atoms in bond
        :return int: GROMACS potential type id
        """
        tipe = cls.gromacs_type_ids[natoms - 2]
        if tipe is None:
            raise TypeError(
                f"The functional form {cls.__name__} does not have a defined GROMACS "
                f"potential type when used with {natoms} atoms.")

        return tipe


class Harmonic(FunctionalForm):
    """Simple harmonic potential.

    See http://manual.gromacs.org/documentation/current/reference-manual/functions/bonded-interactions.html#harmonic-potential  # noqa
    """
    # TODO: Consider whether to use improper (type 2) instead, it is actually harmonic
    gromacs_type_ids = (1, 1, 1)

    def fconst(self, values: typing.Iterable[float], temp: float) -> float:
        rt = 8.314 * temp / 1000.  # pylint: disable=invalid-name
        var = self._variance_func(values)
        return rt / var


class CosHarmonic(FunctionalForm):
    """Cosine based angle potential.

    See http://manual.gromacs.org/documentation/current/reference-manual/functions/bonded-interactions.html#cosine-based-angle-potential  # noqa

    Uses the transformation in eqn 20 of the above source.
    """
    gromacs_type_ids = (None, 2, None)

    def fconst(self, values: typing.Iterable[float], temp: float) -> float:
        rt = 8.314 * temp / 1000.  # pylint: disable=invalid-name
        mean = self.eqm(values, temp)
        var = self._variance_func(values)
        return rt / (math.sin(mean)**2 * var)


class MartiniDefaultLength(FunctionalForm):
    """Dummy functional form which returns a fixed force constant."""
    gromacs_type_ids = (1, None, None)

    def fconst(self, values: typing.Iterable[float], temp: float) -> float:
        return 1250.


class MartiniDefaultAngle(FunctionalForm):
    """Dummy functional form which returns a fixed force constant."""
    gromacs_type_ids = (None, 2, None)

    def fconst(self, values: typing.Iterable[float], temp: float) -> float:
        return 25.


class MartiniDefaultDihedral(FunctionalForm):
    """Dummy functional form which returns a fixed force constant."""
    gromacs_type_ids = (None, None, 1)

    def fconst(self, values: typing.Iterable[float], temp: float) -> float:
        return 50.
