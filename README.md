[![DOI](https://zenodo.org/badge/23794148.svg)](https://zenodo.org/badge/latestdoi/23794148)
[![Travis CI build status](https://travis-ci.org/PCMSolver/pcmsolver.svg?branch=release%2F1.2.Z)](https://travis-ci.org/PCMSolver/pcmsolver)
[![Documentation Status](https://readthedocs.org/projects/pcmsolver/badge/?version=stable)](http://pcmsolver.readthedocs.org/en/latest/?badge=latest)
[![Coverage Status](https://codecov.io/gh/PCMSolver/pcmsolver/branch/release%2F1.2.Z/graph/badge.svg)](https://codecov.io/gh/PCMSolver/pcmsolver)

PCMSolver
=========

An API for the Polarizable Continuum Model. Copyright [Roberto Di Remigio](mailto:roberto.d.remigio@uit.no),
[Luca Frediani](mailto:luca.frediani@uit.no) and [contributors](https://github.com/PCMSolver/pcmsolver/blob/release/1.2.Z/AUTHORS.md)

- [Project website](https://github.com/PCMSolver/pcmsolver)
- [Changelog](CHANGELOG.md)
- [Documentation](http://pcmsolver.readthedocs.io)
- [Build and test history](https://travis-ci.org/PCMSolver/pcmsolver/builds)
- Licensed under [LGPLv3](LICENSE)
- CMake infrastructure managed *via* [Autocmake](http://autocmake.readthedocs.io/)

Continuous integration builds
=============================

All CI builds are triggered by push events to any branch.
Travis CI runs release builds using [ccache](https://ccache.samba.org/) to speed up compilation.

Ubuntu 14.04 LTS 64-bit with CMake 3.5.1 and Boost 1.54.0 this is the
environment offered by [Travis CI](https://travis-ci.org).
Python and Python packages are installed and managed _via_ [Pipenv](http://pipenv.readthedocs.io/en/latest/),
using the `Pipfile` and `Pipfile.lock` files. The following
compilers are used:

1. GCC 4.8, Python 2.7 This build generates _both_ the shared and static
   libraries, linking executables to the former.
2. GCC 6.3.0, Python 2.7 This build generates _only_ the static library.
3. Clang 3.5, GFortran 4.8, Python 3.5 This build generates _both_ the shared and static
   libraries, linking executables to the former.
4. GCC 4.8, Python 2.7 This is a _debug_ build generating _both_ the shared and static
   libraries, linking executables to the former. The build is run with
   coverage analysis for submission to [Codecov](https://codecov.io).
