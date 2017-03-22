[![Build Status](https://travis-ci.org/jag1g13/pycgtool.svg?branch=master)](https://travis-ci.org/jag1g13/pycgtool) [![Documentation Status](https://readthedocs.org/projects/pycgtool/badge/?version=master)](http://pycgtool.readthedocs.io/en/master/?badge=master) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.259330.svg)](https://doi.org/10.5281/zenodo.259330)


# PyCGTOOL
Please see http://pycgtool.readthedocs.io/en/master/ for full documentation.


A Python program for automated generation of coarse-grained molecular dynamics models from atomistic simulation trajectories.

The aim of this project is to provide a tool to aid in parametrising coarse-grained (CG) molecular mechanics models.  PyCGTOOL generates coarse-grained models from atomistic simulation trajectories using a user-provided mapping.  Equilibrium values and force constants of bonded terms are calculated by Boltzmann Inversion of bond distributions collected from the input trajectory.

Alternatively map-only mode (behaving similarly to MARTINIZE) may be used to generate initial coordinates to use with existing CG topologies such as the MARTINI lipid models.  For instance, a pre-equilibrated atomistic membrane may be used to create starting coordinates for a MARTINI membrane simulation.

PyCGTOOL makes it easy to test multiple variations in mapping and bond topology by making simple changes to the config files.

This version has several advantages over the original C++ implementation CGTOOL:
* PyCGTOOL is able to run anywhere the necessary library dependencies are available (all available from pip)
* Does not require that residues are present in contiguous sorted blocks
* May map multiple residues with a single pass
* Support for polymers such as DNA or proteins making use of GROMACS' pdb2gmx
* Much more automated testing ensures that regressions will be identified quickly

If you experience problems or wish to see a new feature added please [file an issue](https://github.com/jag1g13/pycgtool/issues).

If you find this useful, please cite as : `Graham, J. (2017). PyCGTOOL, https://doi.org/10.5281/zenodo.259330`

## Usage
Input to PyCGTOOL is GROMACS GRO and XTC files, along with two custom files: MAP and BND.  These files provide the atomistic-to-CG mapping and bonded topology respectively.  Example files are present in the [test/data](https://github.com/jag1g13/pycgtool/tree/master/test/data) directory.  The format of these files is described in the [full documentation](http://pycgtool.readthedocs.io/en/master/).

To run PyCGTOOL:
`pycgtool.py -g <GRO file> -x <XTC file> -m <MAP file> -b <BND file>`

To run PyCGTOOL in map-only mode:
`pycgtool.py -g <GRO file> -m <MAP file>`

To see the help text:
`pycgtool.py -h`

For more information, see [the tutorial](https://pycgtool.readthedocs.io/en/master/tutorial.html).
It is important to perform validation of any new parameter set; a brief example is present at the end of the tutorial.

## Requirements
PyCGTOOL requires:
* Python3
* [NumPy](http://www.numpy.org/)
* [simpletraj](https://github.com/arose/simpletraj)

The bundled test code may be run using your preferred Python testing frontend although py.test or nose2 is recommended.
All library dependencies may be installed from pip using the command `pip install -r requirements.txt`

This program is a reimplementation of the earlier [CGTOOL](https://bitbucket.org/jag1g13/cgtool).
