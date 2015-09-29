[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.11910.png)](http://dx.doi.org/10.5281/zenodo.11910)
[![Build Status](https://travis-ci.org/PCMSolver/pcmsolver.svg?branch=release)](https://travis-ci.org/PCMSolver/pcmsolver)
[![Magnum CI build status](https://magnum-ci.com/status/9207aa29405095b0b7aef0cd809ed6c2.png?branch=master)](https://magnum-ci.com/builds)
[![Coverage Status](https://coveralls.io/repos/PCMSolver/pcmsolver/badge.svg?branch=release)](https://coveralls.io/r/PCMSolver/pcmsolver?branch=release)
[![Coverity Scan Build](https://scan.coverity.com/projects/3046/badge.svg)](https://scan.coverity.com/projects/3046)

PCMSolver
=========

An API for the Polarizable Continuum Model.

- [Project website](https://gitlab.com/PCMSolver/pcmsolver)
- [Documentation](http://pcmsolver.github.io/pcmsolver-doc)
- [Build and test history](https://travis-ci.org/PCMSolver/pcmsolver/builds)
- [Nightly build dashboard](https://testboard.org/cdash/index.php?project=PCMSolver)
- Version 1.0.3 available
- Licensed under [LGPLv3](../release/COPYING.LESSER)
- CMake infrastructure managed *via* [Autocmake](http://autocmake.readthedocs.org/)

Primary test environments
=========================

All builds force custom build of the needed Boost libraries, except when
stated otherwise.

Continuous integration builds
-----------------------------

The Magnum CI builds are run on push events to any branch, while those
on Travis CI only when pushing to the master/release branch.

- Ubuntu 12.04 LTS 64-bit. GCC 4.6.3, Python 2.7.3, CMake 2.8.7
  This is the environment offered by [Magnum CI](https://magnum-ci.com)
- Ubuntu 12.04 LTS 64-bit. GCC 4.6.3, Python 2.7.3, CMake 2.8.7
  This is the environment offered by [Travis CI](https://travis-ci.org)
- Ubuntu 12.04 LTS 64-bit. Clang 3.4, GFortran 4.6.3, Python 2.7.3, CMake 2.8.7
  This is the environment offered by [Travis CI](https://travis-ci.org)

Nightly builds
--------------

- CentOS 6.6. Intel 12.1.2, Python 2.7.3, CMake 3.1.0
- CentOS 6.6. Intel 13.0, Python 2.7.3, CMake 3.1.0
- CentOS 6.6. Intel 13.4, Python 2.7.3, CMake 3.1.0
- CentOS 6.6. Intel 14.0, Python 2.7.3, CMake 3.1.0
- CentOS 6.6. GCC 4.7.2, Python 2.7.3, CMake 3.1.0
- OS X 10.10.5 Yosemite. LLVM 7.0.0 , GFortan 5.2.0, Python 2.7.10, CMake 3.3.2
  Uses Boost 1.58.0 from the Homebrew repositories.
