import abc
import math

from pycgtool.util import SimpleEnum


class FunctionalForms(object):
    FormsEnum = SimpleEnum.enum("FormsEnum")

    @classmethod
    def refresh(cls):
        enum_dict = cls.FormsEnum.as_dict()
        for subclass in FunctionalForm.__subclasses__():
            name = subclass.__name__
            if name not in cls.FormsEnum:
                enum_dict[name] = subclass

        cls.FormsEnum = SimpleEnum.enum_from_dict("FormsEnum", enum_dict)

    def __init__(self, temp):
        self._temp = temp

    def __getattr__(self, item):
        return type(self).FormsEnum[item].value.boltzmann_inversion

    def __repr__(self):
        return "<FunctionalForms: {0} defined>".format(len(self))

    def __len__(self):
        return len(type(self).FormsEnum)

    def __contains__(self, item):
        return item in type(self).FormsEnum


class FunctionalForm(object, metaclass=abc.ABCMeta):
    @staticmethod
    @abc.abstractstaticmethod
    def boltzmann_inversion(mean, var, temp):
        pass


class Harmonic(FunctionalForm):
    @staticmethod
    def boltzmann_inversion(mean, var, temp):
        rt = 8.314 * temp / 1000.
        return rt / var


class CosHarmonic(FunctionalForm):
    @staticmethod
    def boltzmann_inversion(mean, var, temp):
        rt = 8.314 * temp / 1000.
        rad2 = math.pi * math.pi / (180. * 180.)
        return rt / (math.sin(math.radians(mean))**2 * var * rad2)


class MartiniDefaultLength(FunctionalForm):
    @staticmethod
    def boltzmann_inversion(mean, var, temp):
        return 1250.


class MartiniDefaultAngle(FunctionalForm):
    @staticmethod
    def boltzmann_inversion(mean, var, temp):
        return 25.

FunctionalForms.refresh()


