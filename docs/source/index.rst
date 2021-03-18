.. PyCGTOOL documentation master file, created by
   sphinx-quickstart on Tue Feb 23 21:47:28 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to PyCGTOOL's documentation!
====================================

.. toctree::
   :maxdepth: 2

   self
   tutorial
   mapping-only
   autoapi/index


Features
--------
PyCGTOOL provides a means to quickly and easily generate coarse-grained molecular dynamics models within the MARTINI framework from all-atom or united-atom simulation trajectories.

A user defined mapping is applied to the input trajectory and bonded terms (lengths, angles and dihedrals) are measured.  From these measurements, equilibrium values and force constants are found and a GROMACS topology is created for the target molecules.

Requirements
------------
PyCGTOOL requires:

- Python 3.6 or greater
- Numpy (http://www.numpy.org/)
- MDTraj (http://mdtraj.org/1.7.2/)
- Scipy (https://www.scipy.org/)
- Cython (http://cython.org/)

Optional:

- Python testing framework (e.g. py.test)
- Sphinx to generate documentation yourself (http://www.sphinx-doc.org/en/stable/)

Basic Usage
-----------
PyCGTOOL requires four input files to generate a coarse-grained model:

- AA simulation topology file (e.g. PDB, GRO, etc.)
- AA simulation trajectory file (e.g. XTC, DCD, etc.)
- PyCGTOOL mapping definition file
- PyCGTOOL bond definition file

The program is called by::

    pycgtool <topology file> <trajectory file> -m <MAP file> -b <BND file>

Example mapping and bond definition files are present in the ``test/data`` directory.  Their format is explained below.

After running PyCGTOOL two files, ``out.gro`` and ``out.itp`` will be created.  The gro file contains the mapped coarse-grain coordinates with every molecule for which a mapping was provided.  The itp file contains the parameters for each molecule type.

It is important to perform validation of any new parameter set.
This is typically done by comparing properties between the reference simulation and simulations using the new CG model.
In the tutorial we use the radius of gyration, but there are many other suitable properties, depending on the class of molecule being parametrised.

Mapping / Bond Definition Files
-------------------------------
The mapping and bond definition input files use a format similar to the GROMACS itp/top format.

Mapping Definition
~~~~~~~~~~~~~~~~~~
An example of mapping definition file for the monosaccharide allose taken from ``test/data/sugar.map`` is shown below.

Molecule names (as present in the gro coordinate file) are used as section headers inside square brackets.  Each of the following lines describes a single coarse-grained bead mapping.  The items on a line are: the name of the bead, its MARTINI bead type, optionally the bead charge, and a list of all the atom names it should contain.  All items on a line are whitespace separated.  Multiple molecules may be specified in their own sections.
It is not recommended to provide a mapping for water since MARTINI water combines four molecules into a single bead which is not yet supported by PyCGTOOL.
Note that bead charges in the MARTINI framework are by convention integers and are used only for formally charged functional groups.  An example of a molecule mapping using charges can be found in ``test/data/dppc.map``. ::

   ; comments begin with a semicolon
   [ALLA]
   C1 P3 C1 O1
   C2 P3 C2 O2
   C3 P3 C3 O3
   C4 P3 C4 O4
   C5 P2 C5 C6 O6
   O5 P4 O5

   [SOL]
   W P4 OW HW1 HW2

The Martini 3.0 force field has now introduced the use of virtual sites for polycylic compounds or polysaccharides
in order to improve numerical stability for highly constrained structures. Below is an example of how virtual sites
can be included in the mapping file for naphthalene (see ``/test/data/martini3/naphthalene.map`` ). Similar to the previously
described mapping syntax, a virtual site is defined with a prefix of ``@v`` as follows ``@v [name] [type] [charge] [constructing beads]``.
The constructing beads refer to a list of space delimited coarse grained bead names from which the position of the
virtual site will be calculated. Currently virtual sites can be constructed from either the center of geometry or mass of
the constructing sites via the ``--virtual_map_center`` flag.::

   [ NAPH ]
   R1 TC4 C8 H8 C9 H9
   R2 TC4 C6 H6 C7 H7
   @v R3 TC4 R1 R2 R4 R5
   R4 TC4 C1 H1 C2 H2
   R5 TC4 C3 H3 C4 H4

Bond Definition
~~~~~~~~~~~~~~~
An example bond definition file for the monosaccharide allose taken from ``test/data/sugar.bnd`` is shown below.

As in the mapping definition file, molecule names are used as section headers inside square brackets.  The following lines define bonds lengths, angles and dihedrals between coarse-grained beads.  Each line is a list of bead names, using the names defined in the mapping file.  Two bead names on a line defines a bond length, three defines an angle, and four defines a dihedral.

If no angles are defined for a molecule, PyCGTOOL will construct all angles from the list of bonds.  This may also be enabled for dihedrals via the ``--generate_dihedrals`` flag, but is not recommended as in most cases coarse-grained models do not require dihedrals.  Additionally, any angles inside a triangle of bond lengths are excluded from the output as they often cause simulation stability issues when used in conjunction with LINCS. ::

   ; comments begin with a semicolon
   [ALLA]
   C1 C2
   C2 C3
   C3 C4
   C4 C5
   C5 O5
   O5 C1

   C1 C2 C3
   C2 C3 C4
   C3 C4 C5
   C4 C5 O5
   C5 O5 C1
   O5 C1 C2

   C1 C2 C3 C4
   C2 C3 C4 C5
   C3 C4 C5 O5
   C4 C5 O5 C1
   C5 O5 C1 C2
   O5 C1 C2 C3

Advanced Usage
--------------

Modes
~~~~~

PyCGTOOL performs several other functions which may be useful in the testing and use of coarse-grained models.

Mapping Only
............

Mapping-only mode simply converts an input atomistic coordinate file into its coarse-grained representation.
For full detail see: :doc:`mapping-only`.

This mode may be invoked by::

   pycgtool <topology file> -m <MAP file>


Measure Only
............

Measure-only mode may be used to aid in the testing of a coarse-grained model by making measurements of bonds from a true coarse-grained simulation trajectory.
These bond measurements may be compared directly to those collected from the pseudo-coarse-grained trajectory used to generate the model.
This mode may be invoked by::

   pycgtool <topology file> <trajectory file> -b <BND file>


Advanced Options
~~~~~~~~~~~~~~~~

==================   ==========================================   =======================
Option               Description                                  Values
==================   ==========================================   =======================
output_name          Base name of output files                    **out**, any string
output               Coordinate output format                     **gro**
output_xtc           Should a pseudo-CG XTC be created            **False**, True
map_center           Mapping method                               **geom**, mass
virtual_map_center   Virtual site mapping method                  **geom**, mass
constr_threshold     Convert stiff bonds to constraints over      **100000**, any number
dump_measurements    Whether to output bond measurements          **False**, True
dump_n_values        How many measurements to output              **10000**, any number
output_forcefield    Output a GROMACS forcefield directory?       **False**, True
temperature          Temperature of reference simulation          **310**, any number
default_fc           Use default MARTINI force constants?         **False**, True
generate_angles      Generate angles from bonds                   **True**, False
generate_dihedrals   Generate dihedrals from bonds                **False**, True
==================   ==========================================   =======================
