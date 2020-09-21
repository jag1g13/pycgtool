Mapping Only Mode
=================

PyCGTOOL may be used in mapping only mode to simply convert an all-atom (AA) coordinate file into its coarse-grained (CG) representation.
In this respect it functions similarly to the MARTINI tool `Martinize`.
The main uses of the functionality are:

- Enable direct comparison between AA and CG simulations
- Use a pre-equilibrated AA coordinate file as the starting point for a CG simulation

In order to perform the AA to CG mapping, a PyCGTOOL mapping file is required.
For guidance on the creation of this file see the full :doc:`tutorial`.

To perform a mapping from a single snapshot `.gro` file use the command::

    pycgtool.py -g <GRO file> -m <MAP file>

To perform a mapping of a complete `.xtc` trajectory use the command::

    pycgtool.py -g <GRO file> -x <XTC file> -m <MAP file>

But note that mapping an entire trajectory requires that the optional dependency `MDTraj` is installed.
