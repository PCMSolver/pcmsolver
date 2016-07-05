[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.11910.png)](http://dx.doi.org/10.5281/zenodo.11910)
[![Travis CI Build Status](https://travis-ci.org/PCMSolver/pcmsolver.svg?branch=master)](https://travis-ci.org/PCMSolver/pcmsolver)
[![Magnum CI build status](https://magnum-ci.com/status/9207aa29405095b0b7aef0cd809ed6c2.png?branch=master)](https://magnum-ci.com/builds)
[![AppVeyor CI Build status](https://ci.appveyor.com/api/projects/status/qolwog8ql2gefebr?svg=true)](https://ci.appveyor.com/project/robertodr/pcmsolver)
[![Documentation Status](https://readthedocs.org/projects/pcmsolver/badge/?version=latest)](http://pcmsolver.readthedocs.org/en/latest/?badge=latest)
[![Coverage Status](https://coveralls.io/repos/PCMSolver/pcmsolver/badge.svg?branch=release)](https://coveralls.io/r/PCMSolver/pcmsolver?branch=release)
[![Coverity Scan Build](https://scan.coverity.com/projects/3046/badge.svg)](https://scan.coverity.com/projects/3046)

PCMSolver
=========

An API for the Polarizable Continuum Model.

- [Project website](https://gitlab.com/PCMSolver/pcmsolver)
- [Changelog](../master/CHANGELOG.md)
- [Documentation](http://pcmsolver.readthedocs.org)
- [Build and test history](https://travis-ci.org/PCMSolver/pcmsolver/builds)
- [Nightly build dashboard](https://testboard.org/cdash/index.php?project=PCMSolver)
- Version 1.1.4 available
- Licensed under [LGPLv3](../master/COPYING.LESSER)
- CMake infrastructure managed *via* [Autocmake](http://autocmake.readthedocs.org/)

Primary test environments
=========================

All builds force custom build of the needed Boost libraries, except when
stated otherwise.

Continuous integration builds
-----------------------------

The Magnum CI builds are run on push events to any branch, while those
on Travis CI only when pushing to the `master` branch.
All Travis CI builds on master use ccache to speed up execution.

- Ubuntu 12.04 LTS 64-bit. GCC 4.6, Python 2.7.3, CMake 3.4.2
  This is the environment offered by [Magnum CI](https://magnum-ci.com)
- Ubuntu 12.04 LTS 64-bit with Python 2.7.3, CMake 3.3.2 and Boost 1.55.0
  this is the environment offered by [Travis CI](https://travis-ci.org) pulling
  in various PPA. The following compilers are used, both in release and debug:

  1. GCC 4.6
  2. GCC 4.7
  3. GCC 4.8
  4. GCC 4.9
  5. GCC 5.1, with and without coverage analysis in debug mode
  6. Clang 3.5 and GFortran 4.6
  7. Clang 3.6 and GFortran 4.6
  8. Clang 3.7 and GFortran 4.6
  9. Clang 3.8 and GFortran 4.6

- Mac OS X 10.9.5 with Python 2.7.10, CMake 3.2.3 and Boost 1.58.0
  this is the environment offered by [Travis CI](https://travis-ci.org)
  The following compilers are used, both in release and debug:

  1. XCode 6.4 with Clang and GFortran 5.2
  2. XCode 6.4 with GCC 5.2
  3. XCode 7.0 with Clang and GFortran 5.2
  4. XCode 7.0 with GCC 5.2

The build needed for submission to [Coverity scan](https://scan.coverity.com/)
is triggered by pushes to the `coverity_scan` branch. It is run on
Ubuntu 12.04 LTS 64-bit with Python 2.7.3, CMake 3.3.2 and Boost 1.55.0
this is the environment offered by [Travis CI](https://travis-ci.org) pulling
in various PPA. GCC 5.1 is used, in debug mode.

Nightly builds
--------------

- CentOS 6.6. Intel 12.1.2, Python 2.7.3, CMake 3.1.0
- CentOS 6.6. Intel 13.0, Python 2.7.3, CMake 3.1.0
- CentOS 6.6. Intel 13.4, Python 2.7.3, CMake 3.1.0
- CentOS 6.6. Intel 14.0, Python 2.7.3, CMake 3.1.0
- CentOS 6.6. Intel 15.0, Python 2.7.9, CMake 3.2.2.
  Uses Boost 1.58.0
- CentOS 6.6. GCC 4.4.7, Python 2.7.3, CMake 3.1.0
- CentOS 6.6. GCC 4.7.2, Python 2.7.3, CMake 3.1.0
- CentOS 6.6. GCC 4.9.1, Python 2.7.3, CMake 3.1.0
- OS X 10.10.5 Yosemite. LLVM 7.0.0 , GFortran 5.2.0, Python 2.7.10, CMake 3.3.2
  Uses Boost 1.58.0 from the Homebrew repositories.
