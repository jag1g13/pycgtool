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

    @staticmethod
    @abc.abstractstaticmethod
    def fconst(values, temp):
        """
        Calculate force constant.
        Abstract static method to be defined by all functional forms.

        :param values: Measured internal coordinate values from which to calculate force constant
        :param temp: Temperature of simulation
        :return: Calculated force constant
        """
        pass


class Harmonic(FunctionalForm):
    @staticmethod
    def fconst(values, temp):
        rt = 8.314 * temp / 1000.
        var = np.nanvar(values)
        return rt / var


class CosHarmonic(FunctionalForm):
    @staticmethod
    def fconst(values, temp):
        rt = 8.314 * temp / 1000.
        mean = CosHarmonic.eqm(values, temp)
        var = np.nanvar(values)
        return rt / (math.sin(mean)**2 * var)


class MartiniDefaultLength(FunctionalForm):
    @staticmethod
    def fconst(values, temp):
        return 1250.


class MartiniDefaultAngle(FunctionalForm):
    @staticmethod
    def fconst(values, temp):
        return 25.


class MartiniDefaultDihedral(FunctionalForm):
    @staticmethod
    def fconst(values, temp):
        return 50.
