import numpy as np

from .util import stat_moments
from .parsers.cfg import CFG


class Bond:
    __slots__ = ["atoms", "atom_numbers", "values", "eqm", "fconst"]

    def __init__(self, atoms=None, atom_numbers=None):
        self.atoms = atoms
        self.atom_numbers = atom_numbers
        self.values = []

    def boltzmann_invert(self, temp=310):
        # TODO needs tests
        mean, var = stat_moments(self.values)

        rt = 8.314 * temp / 1000.
        rad2 = np.pi * np.pi / (180. * 180.)
        conv = {2: lambda: rt / var,
                3: lambda: rt / (np.sin(np.radians(mean))**2 * var * rad2),
                4: lambda: rt / (var * rad2)}

        self.eqm = mean
        try:
            self.fconst = conv[len(self.atoms)]()
        except FloatingPointError:
            self.fconst = 0

    def __repr__(self):
        try:
            return "<Bond containing atoms {0} with r_0 {1:.3f} and force constant {2:.3e}>".format(
                ", ".join(self.atoms), self.eqm, self.fconst
            )
        except AttributeError:
            return "<Bond containing atoms {0}>".format(", ".join(self.atoms))


def angle(a, b, c=None):
    if c is None:
        c = np.cross(a, b)
    dot = np.dot(a, b)
    abscross = np.sqrt(np.dot(c, c))
    return np.arctan2(abscross, dot)


class BondSet:
    def __init__(self, filename=None):
        self._molecules = {}

        if filename is not None:
            with CFG(filename) as cfg:
                for mol in cfg:
                    self._molecules[mol.name] = []
                    molbnds = self._molecules[mol.name]
                    for atomlist in mol:
                        molbnds.append(Bond(atoms=atomlist))

    def _populate_atom_numbers(self, mapping):
        for mol in self._molecules:
            molmap = mapping[mol]
            for bond in self._molecules[mol]:
                ind = lambda name: [bead.name for bead in molmap].index(name)
                bond.atom_numbers = [ind(atom) for atom in bond.atoms]

    def write_itp(self, filename, mapping):
        self._populate_atom_numbers(mapping)

        with open(filename, "w") as itp:
            header = ("; \n"
                      "; Topology prepared automatically using PyCGTOOL \n"
                      "; James Graham <J.A.Graham@soton.ac.uk> 2016 \n"
                      "; University of Southampton \n"
                      "; https://github.com/jag1g13/pycgtool \n"
                      ";")
            print(header, file=itp)

            # Print molecule
            for mol in self._molecules:
                print("\n[ moleculetype ]", file=itp)
                print("{0:4s} {1:4d}".format(mol, 1), file=itp)

                print("\n[ atoms ]", file=itp)
                for i, bead in enumerate(mapping[mol]):
                    print("{0:4d} {1:4s} {2:4d} {3:4s} {4:4s} {5:4d} {6:8.3f}".format(
                        i+1, bead.type, 1, mol, bead.name, 1, bead.charge
                    ), file=itp)

                bonds = [bond for bond in self._molecules[mol] if len(bond.atoms) == 2]
                if len(bonds):
                    print("\n[ bonds ]", file=itp)
                for bond in bonds:
                    print("{0:4d} {1:4d} {2:4d} {3:12.5f} {4:12.5f}".format(
                        bond.atom_numbers[0]+1, bond.atom_numbers[1]+1,
                        1, bond.eqm, bond.fconst
                    ), file=itp)

                bonds = [bond for bond in self._molecules[mol] if len(bond.atoms) == 3]
                if len(bonds):
                    print("\n[ angles ]", file=itp)
                for bond in bonds:
                    print("{0:4d} {1:4d} {2:4d} {3:4d} {4:12.5f} {5:12.5f}".format(
                        bond.atom_numbers[0]+1, bond.atom_numbers[1]+1, bond.atom_numbers[2]+1,
                        1, bond.eqm, bond.fconst
                    ), file=itp)

                bonds = [bond for bond in self._molecules[mol] if len(bond.atoms) == 4]
                if len(bonds):
                    print("\n[ dihedrals ]", file=itp)
                for bond in bonds:
                    print("{0:4d} {1:4d} {2:4d} {3:4d} {4:4d} {5:12.5f} {6:12.5f} {7:4d}".format(
                        bond.atom_numbers[0]+1, bond.atom_numbers[1]+1,
                        bond.atom_numbers[2]+1, bond.atom_numbers[3]+1,
                        1, bond.eqm, bond.fconst, 1
                    ), file=itp)

    def apply(self, frame):
        def calc_length(res, atoms):
            vec = res[atoms[1]].coords - res[atoms[0]].coords
            return np.sqrt(np.sum(vec*vec))

        def calc_angle(res, atoms):
            # for atom in atoms:
            # print(res[atom].coords)
            veca = res[atoms[1]].coords - res[atoms[0]].coords
            vecb = res[atoms[2]].coords - res[atoms[1]].coords
            ret = np.degrees(np.pi - angle(veca, vecb))
            # print(ret)
            return ret

        def calc_dihedral(res, atoms):
            vec1 = res[atoms[1]].coords - res[atoms[0]].coords
            vec2 = res[atoms[2]].coords - res[atoms[1]].coords
            vec3 = res[atoms[3]].coords - res[atoms[2]].coords

            c1 = np.cross(vec1, vec2)
            c2 = np.cross(vec2, vec3)
            c3 = np.cross(c1, c2)

            ang = np.degrees(angle(c1, c2))
            direction = np.dot(vec2, c3)
            return ang if direction > 0 else -ang

        for res in frame:
            mol_meas = self._molecules[res.name]
            for bond in mol_meas:
                try:
                    calc = {2: calc_length,
                            3: calc_angle,
                            4: calc_dihedral}
                    val = calc[len(bond.atoms)](res, bond.atoms)
                    bond.values.append(val)
                except NotImplementedError:
                    pass

    def boltzmann_invert(self):
        for mol in self._molecules:
            for bond in self._molecules[mol]:
                bond.boltzmann_invert()

    def __len__(self):
        return len(self._molecules)

    def __contains__(self, item):
        return item in self._molecules

    def __getitem__(self, item):
        return self._molecules[item]

    def __iter__(self):
        return iter(self._molecules)

