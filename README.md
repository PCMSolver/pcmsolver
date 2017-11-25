[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.11910.png)](http://dx.doi.org/10.5281/zenodo.11910)
[![Travis CI build status](https://travis-ci.org/PCMSolver/pcmsolver.svg?branch=release%2F1.Y)](https://travis-ci.org/PCMSolver/pcmsolver)
[![Documentation Status](https://readthedocs.org/projects/pcmsolver/badge/?version=stable)](http://pcmsolver.readthedocs.org/en/latest/?badge=latest)
[![Coverage Status](https://coveralls.io/repos/PCMSolver/pcmsolver/badge.svg?branch=release%2F1.Y)](https://coveralls.io/r/PCMSolver/pcmsolver?branch=release)
[![Coverity Scan Build](https://scan.coverity.com/projects/3046/badge.svg)](https://scan.coverity.com/projects/3046)

PCMSolver
=========

An API for the Polarizable Continuum Model.

- [Project website](https://github.com/PCMSolver/pcmsolver)
- [Changelog](CHANGELOG.md)
- [Documentation](http://pcmsolver.readthedocs.io)
- [Build and test history](https://travis-ci.org/PCMSolver/pcmsolver/builds)
- [Build dashboard](https://testboard.org/cdash/index.php?project=PCMSolver)
- Version 1.1.11 available
- Licensed under [LGPLv3](LICENSE)
- CMake infrastructure managed *via* [Autocmake](http://autocmake.readthedocs.io/)

Primary test environments
=========================

All builds force custom build of the needed Boost libraries, except when
stated otherwise.

Continuous integration builds
-----------------------------

All CI builds are triggered by push events to any branch.
Travis CI runs release builds using [ccache](https://ccache.samba.org/) to speed up compilation.

- Ubuntu 14.04 LTS 64-bit with CMake 3.5.1 and Boost 1.54.0 this is the
  environment offered by [Travis CI](https://travis-ci.org).
  Python and Python packages are installed and managed _via_ a virtual environment.
  Python dependencies are listed in the `requirements.txt` file. The following
  compilers are used:

  1. GCC 4.6, Python 2.7 This build generates _both_ the shared and static
     libraries, linking executables to the former. The build is run with and
     without coverage analysis, the latter being a _debug_ build.
  2. Clang 3.5, GFortran 4.6, Python 3.5 This build generates _only_ the static
     library.

- Mac OS X 10.11 with CMake 3.6.2 and Boost 1.61.0
  this is the environment offered by [Travis CI](https://travis-ci.org)
  with their Xcode 7.3.1 image.
  The following compilers are used:

  1. Apple Clang 7.3.0, GFortran 5.5.0, Python 2.7 This build generates _only_
     the static library.
  2. GCC 6.4.0, Python 3.6 This build generates _both_ the shared and static
     libraries, linking executables to the former.

The build needed for submission to [Coverity scan](https://scan.coverity.com/)
is triggered by pushes to the `coverity_scan` branch. It is run on
Ubuntu 12.04 LTS 64-bit with Python 2.7, CMake 3.3.2 and Boost 1.55.0
this is the environment offered by [Travis CI](https://travis-ci.org) pulling
in various PPA. GCC 5.1 is used, in debug mode.
