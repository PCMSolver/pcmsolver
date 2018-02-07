[![DOI](https://zenodo.org/badge/23794148.svg)](https://zenodo.org/badge/latestdoi/23794148)
[![Travis CI build status](https://travis-ci.org/PCMSolver/pcmsolver.svg?branch=release%2F1.Y)](https://travis-ci.org/PCMSolver/pcmsolver)
[![Documentation Status](https://readthedocs.org/projects/pcmsolver/badge/?version=stable)](http://pcmsolver.readthedocs.org/en/latest/?badge=latest)
[![Coverage Status](https://codecov.io/gh/PCMSolver/pcmsolver/branch/release%2F1.Y/graph/badge.svg)](https://codecov.io/gh/PCMSolver/pcmsolver)
[![Coverity Scan Build](https://scan.coverity.com/projects/3046/badge.svg)](https://scan.coverity.com/projects/3046)

PCMSolver
=========

An API for the Polarizable Continuum Model. Copyright [Roberto Di Remigio](mailto:roberto.d.remigio@uit.no),
[Luca Frediani](mailto:luca.frediani@uit.no) and [contributors](https://github.com/PCMSolver/pcmsolver/blob/release/1.Y/AUTHORS.md)

- [Project website](https://github.com/PCMSolver/pcmsolver)
- [Changelog](CHANGELOG.md)
- [Documentation](http://pcmsolver.readthedocs.io)
- [Build and test history](https://travis-ci.org/PCMSolver/pcmsolver/builds)
- Version 1.1.12 available
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

1. GCC 4.6, Python 2.7 This build generates _both_ the shared and static
   libraries, linking executables to the former.
2. GCC 6.3.0, Python 2.7 This build generates _only_ the static library.
3. Clang 3.5, GFortran 4.6, Python 3.5 This build generates _both_ the shared and static
   libraries, linking executables to the former.
4. GCC 4.8, Python 2.7 This is a _debug_ build generating _both_ the shared and static
   libraries, linking executables to the former. The build is run with
   coverage analysis for submission to [Codecov](https://codecov.io).

The build needed for submission to [Coverity scan](https://scan.coverity.com/)
is triggered by pushes to the `coverity_scan` branch. It is run on
Ubuntu 12.04 LTS 64-bit with Python 2.7, CMake 3.3.2 and Boost 1.55.0
this is the environment offered by [Travis CI](https://travis-ci.org) pulling
in various PPA. GCC 5.1 is used, in debug mode.
