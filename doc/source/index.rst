.. PyCGTOOL documentation master file, created by
   sphinx-quickstart on Tue Feb 23 21:47:28 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to PyCGTOOL's documentation!
====================================

.. toctree::
   :hidden:

   self
   tutorial
   mapping-only
   Module Documentation <modules>


.. toctree::
   :maxdepth: 2

:doc:`tutorial`

Features
--------
PyCGTOOL provides a means to quickly and easily generate coarse-grained molecular dynamics models within the MARTINI framework from all-atom or united-atom simulation trajectories.

A user defined mapping is applied to the input trajectory and bonded terms (lengths, angles and dihedrals) are measured.  From these measurements, equilibrium values and force constants are found and a GROMACS topology is created for the target molecules.

Requirements
------------
PyCGTOOL requires:

- Python 3
- Numpy (http://www.numpy.org/)
- simpletraj (https://github.com/arose/simpletraj)

Optional:

- MDTraj for pseudo-CG XTC output (http://mdtraj.org/1.7.2/) with own dependencies:

  - Scipy (https://www.scipy.org/)
  - Cython (http://cython.org/)

- Python testing framework (e.g. Nose2, py.test)
- Numba for increased performance (http://numba.pydata.org/)
- Sphinx to generate documentation yourself (http://www.sphinx-doc.org/en/stable/)

Basic Usage
-----------
PyCGTOOL requires four input files to generate a coarse-grained model:

-g       GROMACS gro coordinate file
-x       GROMACS xtc trajectory file
-m       PyCGTOOL mapping definition file
-b       PyCGTOOL bond definition file

The program is called by::

    pycgtool.py -g <GRO file> -x <XTC file> -m <MAP file> -b <BND file>

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

Bond Definition
~~~~~~~~~~~~~~~
An example bond definition file for the monosaccharide allose taken from ``test/data/sugar.bnd`` is shown below.

As in the mapping definition file, molecule names are used as section headers inside square brackets.  The following lines define bonds lengths, angles and dihedrals between coarse-grained beads.  Each line is a list of bead names, using the names defined in the mapping file.  Two bead names on a line defines a bond length, three defines an angle, and four defines a dihedral.

If no angles are defined for a molecule, PyCGTOOL will construct all angles from the list of bonds.  This may also be enabled for dihedrals via the ``--advanced`` flag, but is not recommended as in most cases coarse-grained models do not require dihedrals.  Additionally, any angles inside a triangle of bond lengths are excluded from the output as they often cause simulation stability issues when used in conjunction with LINCS. ::

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

   pycgtool.py -g <GRO file> -m <MAP file>


Measure Only
............

Measure-only mode may be used to aid in the testing of a coarse-grained model by making measurements of bonds from a true coarse-grained simulation trajectory.
These bond measurements may be compared directly to those collected from the pseudo-coarse-grained trajectory used to generate the model.
This mode may be invoked by::

   pycgtool.py -g <GRO file> -x <XTC file> -b <BND file>


Advanced Options
~~~~~~~~~~~~~~~~
By passing the flag ``--advanced`` to PyCGTOOL several advanced options are accessible.  The arrow keys may be used to navigate through the menu.  Enter selects an option to be edited, or if the option is boolean toggles it.  Once you have edited an option press enter again.  When all options are satisfactory, press q to proceed.

==================   ==========================================   =======================
Option               Description                                  Values
==================   ==========================================   =======================
output_name          Base name of output files                    **out**, any string
output               Coordinate output format                     **gro**
output_xtc           Should a pseudo-CG XTC be created            **False**, True
map_only             Run in mapping-only mode                     **False**, True
map_center           Mapping method                               **geom**, mass
constr_threshold     Convert stiff bonds to constraints over      **100000**, any number
dump_measurements    Whether to output bond measurements          **False**, True
dump_n_values        How many measurements to output              **10000**, any number
output_forcefield    Output a GROMACS forcefield directory?       **False**, True
temperature          Temperature of reference simulation          **310**, any number
default_fc           Use default MARTINI force constants?         **False**, True
generate_angles      Generate angles from bonds                   **True**, False
generate_dihedrals   Generate dihedrals from bonds                **False**, True
==================   ==========================================   =======================

Indexes
=======

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

