# PyCGTOOL

[![License](https://img.shields.io/github/license/jag1g13/pycgtool.svg)](LICENSE)
[![Python package](https://github.com/jag1g13/pycgtool/actions/workflows/python-package.yml/badge.svg?branch=dev)](https://github.com/jag1g13/pycgtool/actions)
[![Documentation](https://readthedocs.org/projects/pycgtool/badge/?version=dev)](http://pycgtool.readthedocs.io/en/dev)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.598143.svg)](https://doi.org/10.5281/zenodo.598143)
[![PyPi Version](https://img.shields.io/pypi/v/pycgtool.svg)](https://pypi.python.org/pypi/pycgtool/)
[![Downloads](https://pepy.tech/badge/pycgtool)](https://pepy.tech/project/pycgtool)

Generate coarse-grained molecular dynamics models from atomistic trajectories.

The aim of this project is to provide a tool to aid in parametrising coarse-grained (CG) molecular mechanics models.
PyCGTOOL generates coarse-grained models from atomistic simulation trajectories using a user-provided mapping. 
Equilibrium values and force constants of bonded terms are calculated by Boltzmann Inversion of bond distributions collected from the input trajectory.

Alternatively map-only mode (behaving similarly to MARTINIZE) may be used to generate initial coordinates to use with existing CG topologies such as the MARTINI lipid models.
For instance, a pre-equilibrated atomistic membrane may be used to create starting coordinates for a MARTINI membrane simulation.

PyCGTOOL makes it easy to test multiple variations in mapping and bond topology by making simple changes to the config files.

If you find PyCGTOOL useful, please cite our JCIM paper (https://doi.org/10.1021/acs.jcim.7b00096) and the code itself (https://doi.org/10.5281/zenodo.598143).

```bibtex
@article{Graham2017,
   author = {James A. Graham and Jonathan W. Essex and Syma Khalid},
   doi = {10.1021/acs.jcim.7b00096},
   issn = {1549-9596},
   issue = {4},
   journal = {Journal of Chemical Information and Modeling},
   month = {4},
   pages = {650-656},
   title = {PyCGTOOL: Automated Generation of Coarse-Grained Molecular Dynamics Models from Atomistic Trajectories},
   volume = {57},
   url = {https://pubs.acs.org/doi/10.1021/acs.jcim.7b00096},
   year = {2017},
}
```

## Install

PyCGTOOL requires Python 3.6 or higher and may be installed using `pip`:
```
pip install pycgtool
```

Alternatively, you may download a pre-packaged version for your operating system from the [releases page](https://github.com/jag1g13/pycgtool/releases) on GitHub.
These pre-packaged versions include all dependencies and should be suitable in cases where you cannot install packages with `pip`.
**Warning**: This installation method is not extensively tested - installing via `pip` should be prefered in most cases.

### MDTraj on macOS

On some versions macOS, with some versions of the Clang compiler, MDTraj may fail to load GROMACS XTC simulation trajectories.
If you encounter this issue, make sure you have the latest version of MDTraj.

For more information see [MDTraj/#1572](https://github.com/mdtraj/mdtraj/issues/1572).

## Usage

Input to PyCGTOOL is an atomistic simulation trajectory in the form of a topology (e.g. PDB, GRO, etc.) and a trajectory file (e.g. XTC, DCD, etc.), along with two custom files: MAP and BND.
These files provide the atomistic-to-CG mapping and bonded topology respectively.

Example files are present in the [test/data](https://github.com/jag1g13/pycgtool/tree/master/test/data) directory.
The format of these files is described in the [full documentation](https://pycgtool.readthedocs.io/en/master/index.html).

For more information, see [the tutorial](https://pycgtool.readthedocs.io/en/master/tutorial.html).
It is important to perform validation of any new parameter set - a brief example is present at the end of the tutorial.

For a full list of options, see the [documentation](https://pycgtool.readthedocs.io/en/master/index.html) or use:
```
pycgtool -h
```

### Generate a Model

To generate a CG model from an atomistic simulation:
```
pycgtool <topology file> <trajectory file> -m <MAP file> -b <BND file>
```

### Map Only

To use PyCGTOOL to convert a set of atomistic simulation coordinates to CG coordinates:
```
pycgtool <topology file> -m <MAP file>
```

Or to convert a complete simulation trajectory:
```
pycgtool <topology file> <trajectory file> -m <MAP file>
```

## Maintainers

James Graham ([@jag1g13](https://github.com/jag1g13))

## Contributing

If you experience problems using PyCGTOOL or wish to see a new feature added please [open an issue](https://github.com/jag1g13/pycgtool/issues/new).

To help develop PyCGTOOL, you can create a fork of this repository, clone your fork and install PyCGTOOL in development mode using [Poetry](https://python-poetry.org/):
```
poetry install
```

This will result in an editable mode install (similar to `pip install -e .`) along with all the necessary runtime and development dependencies.
Testing and linting is handled by [Tox](https://tox.readthedocs.io/en/latest/) - use `tox` to run the full test suite and linter as they are configured in the Continuous Integration pipeline.

When you're ready for your work to be merged, please submit a Pull Request.

## License

[GPL-3.0](LICENSE) © James Graham, University of Southampton
