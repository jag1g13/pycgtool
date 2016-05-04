[![Build Status](https://travis-ci.org/jag1g13/pycgtool.svg?branch=master)](https://travis-ci.org/jag1g13/pycgtool)

# PyCGTOOL
Please see http://pycgtool.readthedocs.io/en/master/ for full documentation.

Python reimplementation of [CGTOOL](https://bitbucket.org/jag1g13/cgtool) performing coarse-grain mapping of molecular dynamics trajectories.

The aim of this project is to provide a tool to aid in parametrising coarse-grained (CG) molecular mechanics models.  PyCGTOOL generates coarse-grained models from atomistic trajectories using a user-provided mapping.  Equilibrium values and force constants of bonded terms are calculated by Boltzmann Inversion of histograms collected from the input trajectory allowing good replication of target properties.

Alternatively it may be used in map-only mode (behaving similarly to MARTINIZE) to generate initial coordinates to use with existing CG topologies such as the MARTINI lipid models.

PyCGTOOL makes it easy to test multiple variations in mapping and bond topology by making simple changes to the config file.

This version has several advantages over the original C++ implementation CGTOOL:
* PyCGTOOL is able to run anywhere the necessary library dependencies are available (all available from pip)
* Does not require that residues are present in contiguous sorted blocks
* Much more automated testing ensures that regressions will be identified quickly
* Work-in-progress support for polymers such as DNA or proteins making use of GROMACS' pdb2gmx

If you experience problems or wish to see a new feature added please [file an issue](https://github.com/jag1g13/pycgtool/issues).

## Usage
Input to PyCGTOOL is GROMACS GRO and XTC files, along with two custom files: MAP and BND.  These files provide the atomistic-to-CG mapping and bonded topology respectively.  Example files are present in the [test/data](https://github.com/jag1g13/pycgtool/tree/master/test/data) directory.  The format of these files is described in the [full documentation] (http://pycgtool.readthedocs.io/en/master/).

To run PyCGTOOL:
`pycgtool.py -g <GRO file> -x <XTC file> -m <MAP file> -b <BND file>`

To run PyCGTOOL in map-only mode:
`pycgtool.py -g <GRO file> -m <MAP file>`

To see the help text:
`pycgtool.py -h`

## Requirements
PyCGTOOL requires:
* Python3
* [NumPy](http://www.numpy.org/)
* [simpletraj](https://github.com/arose/simpletraj)

The bundled test code may be run using your preferred Python testing frontend although nose is recommended.
All library dependencies may be installed from pip using the command `pip install -r requirements.txt`
