import abc
import math

import numpy as np

from pycgtool.util import SimpleEnum


class FunctionalForms(object):
    """
    Class holding list of all defined functional forms for Boltzmann Inversion.

    Creating an instance causes the Enum of functional forms to be updated with
    all new subclasses of FunctionalForm.  These may then be accessed by name,
    either as attributes or using square brackets.
    """
    FormsEnum = SimpleEnum.enum("FormsEnum")

    @classmethod
    def _refresh(cls):
        """
        Update the functional forms Enum to include all new subclasses of FunctionalForm.
        """
        enum_dict = cls.FormsEnum.as_dict()
        for subclass in FunctionalForm.__subclasses__():
            name = subclass.__name__
            if name not in cls.FormsEnum:
                enum_dict[name] = subclass()

        cls.FormsEnum = SimpleEnum.enum_from_dict("FormsEnum", enum_dict)

    def __init__(self, **kwargs):
        self._kwargs = kwargs
        type(self)._refresh()

    def __getattr__(self, item):
        return type(self).FormsEnum[item].value

    def __getitem__(self, item):
        return getattr(self, item)

    def __repr__(self):
        return "<FunctionalForms: {0} defined>".format(len(self))

    def __len__(self):
        return len(type(self).FormsEnum)

    def __contains__(self, item):
        return item in type(self).FormsEnum


class FunctionalForm(object, metaclass=abc.ABCMeta):
    """
    Parent class of any functional form used in Boltzmann Inversion to convert variance to a force constant.

    New functional forms must define a static __call__ method.
    """
    @staticmethod
    def eqm(values, temp):
        """
        Calculate equilibrium value.
        May be overridden by functional forms.

        :param values: Measured internal coordinate values from which to calculate equilibrium value
        :param temp: Temperature of simulation
        :return: Calculated equilibrium value
        """
        return np.nanmean(values)

    @abc.abstractstaticmethod
    def fconst(values, temp):
        """
        Calculate force constant.
        Abstract static method to be defined by all functional forms.

        :param values: Measured internal coordinate values from which to calculate force constant
        :param temp: Temperature of simulation
        :return: Calculated force constant
        """
        raise NotImplementedError

    @abc.abstractproperty
    def gromacs_type_ids(self):
        """
        Return tuple of GROMACS potential type ids when used as length, angle, dihedral.
        
        :return tuple[int]: Tuple of GROMACS potential type ids
        """
        raise NotImplementedError

    @classmethod
    def gromacs_type_id_by_natoms(cls, natoms):
        """
        Return the GROMACS potential type id for this functional form when used with natoms.
        
        :param int natoms: 
        :return int: GROMACS potential type id 
        """
        tipe = cls.gromacs_type_ids[natoms - 2]
        if tipe is None:
            raise TypeError("The functional form {0} does not have a defined GROMACS potential type when used with {1} atoms.".format(cls.__name__, natoms))
        return tipe


class Harmonic(FunctionalForm):
    gromacs_type_ids = (1, 1, 1)  # Consider whether to use improper (type 2) instead, it is actually harmonic

    @staticmethod
    def fconst(values, temp):
        rt = 8.314 * temp / 1000.
        var = np.nanvar(values)
        return rt / var


class CosHarmonic(FunctionalForm):
    gromacs_type_ids = (None, 2, None)

    @staticmethod
    def fconst(values, temp):
        rt = 8.314 * temp / 1000.
        mean = CosHarmonic.eqm(values, temp)
        var = np.nanvar(values)
        return rt / (math.sin(mean)**2 * var)


class MartiniDefaultLength(FunctionalForm):
    gromacs_type_ids = (1, None, None)

    @staticmethod
    def fconst(values, temp):
        return 1250.


class MartiniDefaultAngle(FunctionalForm):
    gromacs_type_ids = (None, 2, None)

    @staticmethod
    def fconst(values, temp):
        return 25.


class MartiniDefaultDihedral(FunctionalForm):
    gromacs_type_ids = (None, None, 1)

    @staticmethod
    def fconst(values, temp):
        return 50.
