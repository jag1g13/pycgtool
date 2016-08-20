PyCGTOOL Tutorial
=================

This tutorial follows the complete process of parametrising a new molecule within the MARTINI forcefield, covering aspects of mapping design, model generation and model validation.
PyCGTOOL is used at multiple stages, showing its use in several different situations.

The molecule chosen as a target for this parametrisation is the :math:`\beta_1` antagonist atenolol.


Atomistic Simulation
--------------------
The reference simulation for the parametrisation of atenolol was performed using the GROMOS 54A7 united atom forcefield with a topology from the `ATB database <https://atb.uq.edu.au/molecule.py?molid=23433>`_.
A single molecule of atenolol was solvated and equilibrated, before collecting a 50 ns trajectory.
Currently, PyCGTOOL is limited to trajectories in the GROMACS XTC format, with a single frame GRO for residue information.

Mapping Design
--------------
Designing a suitable mapping from the atomistic to the coarse-grained representation requires some experience and a degree of `chemical intuition`, but the ease with which the mapping may be modified using PyCGTOOL allows the mapping to be iterated much more quickly.

The process of designing a mapping involves splitting the molecule into fragments, each of which contains approximately four heavy atoms.
Start by finding functional groups such as amides or carboxylates, each of which may become a single bead.
Next, replace phenyl rings with a triangle of three `small` type beads, each of which contains two heavy atoms and has reduced Lennard-Jones size and mass, as compared to the normal four-atom beads.
Finally, divide any remaining parts of the molecule into beads of four heavy atoms as required.
The ideal bead will contain four heavy atoms and be nearly spherical, but this is not always possible.
If any atoms remain after clusters of four have been allocated, it may be required to use a mapping of three-to-one for some beads.

After the atoms have been allocated to beads, determine which beads should be bonded by copying the bonds between their component atoms.
It will probably be the case that there is no obviously best mapping and bond topology, which is not at this point a major issue as multiple mappings can be assessed easily.

Once the mapping and bond topology have been determined, they must be put into a format readable by PyCGTOOL.
This format is as described in the introduction to PyCGTOOL.

Model Generation
----------------

Running the CG Simulation
-------------------------

Model Validation
----------------

