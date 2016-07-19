# Change Log

## [Unreleased]

### Added

### Changed

### Deprecated

### Removed

### Fixed

### Security

## [Version 1.1.5] - 2016-07-19

### Added

- A radii set derived from [Allinger's MM3 model](http://dx.doi.org/10.1016/S0166-1280(09)80008-0)
  can now be chosen to build the van der Waals cavity surface.
  Notice that the values reported in the original paper are **divided by** 1.2, to match the
  default radii set used in [ADF](https://www.scm.com/doc/ADF/Input/COSMO.html)
  The closest match to ADF can be obtained by using CPCM as solver, Allinger's radii and setting
  the scaling of radii to false.

## [Version 1.1.4] - 2016-07-05

### Changed

- The `CPCMSolver` object now stores the scaled, Hermitian, symmetry-adapted S matrix.
  Polarization weights are then directly computed from the incoming MEP.
- The `IEFSolver` object now stores the non-Hermitian, symmetry-adapted T and R matrices.
  The explicit calculation of the inverse of T is thus avoided.
  However, two square matrices of size equal to the cavity size are stored instead
  of just one. To obtain the polarization weights _two_ linear systems of equations are solved.
  A [partially pivoted LU decomposition](http://eigen.tuxfamily.org/dox/classEigen_1_1PartialPivLU.html)
  is used to solve the linear system(s).
  The strategy used in v1.1.3 suffered from a reduced numerical accuracy, due to the fact that
  the polarization weights were not correctly defined.

### Removed

- The `hermitivitize` function will only work correctly on matrices. This
  reverts modifications in the previous release.

## [Version 1.1.3] - 2016-07-03

### Changed

- The `PEDRA.OUT` cavity generator log now reports the initial _and_ final
  lists of spheres. The final list only contains those spheres that were
  actually tesselated and, possibly, the added spheres.
- For all solvers, the symmetrization needed to obtain the polarization weights
  happens directly on the computed charges, instead of symmetrizing the system
  matrices.
- The `IEFSolver` object stores the unsymmetrized T^-1R matrices.
  A [partially pivoted LU decomposition](http://eigen.tuxfamily.org/dox/classEigen_1_1PartialPivLU.html)
  is used to compute T^-1R.
  A [robust Cholesky decomposition](http://eigen.tuxfamily.org/dox/classEigen_1_1LDLT.html) is
  used to form the R matrix in the anisotropic IEF case.
- The `CPCMSolver` object now stores the scaled, unsymmetrized S matrix. The
  explicit calculation and storage of its inverse is thus avoided.
  A [robust Cholesky decomposition](http://eigen.tuxfamily.org/dox/classEigen_1_1LDLT.html) is
  used to solve the linear equation system.
- The `hermitivitize` function can now correctly symmetrize vectors.

### Fixed

- A fix for the initialization of the explicit list of spheres when running the
  standalone executable. The bug prevented the generation of the correct
  `Molecule` object, with subsequent failure in the cavity generator.
- A memory leak occuring in the cavity generator PEDRA was fixed. This was uncovered by @shoefener
  and manifested only with considerably large cavities (> 200 input spheres)

### Removed

- The function `CPCMMatrix` in the `SolverImpl.hpp` header file is no longer available.

## [Version 1.1.2] - 2016-05-31

### Fixed

- Signatures for strings in Fortran90 bindings. They have now the proper
  C interoperable type `character(kind=c_char, len=1) :: label(lbl_size)`.
  For the host this means that surface function labels will have to be declared
  as character arrays, for example: `character :: label(7) = (/'T', 'o', 't', 'M', 'E', 'P', char(0)/)`

### Changed

- More informative error messages for runtime crashes caused by access to
  surface functions.
- The signatures for the interface functions now accept and/or return `int` (`c_int`)
  instead of `size_t` (`c_size_t`). This simplifies interfacing with Fortran hosts.

## [Version 1.1.1] - 2016-03-10

### Added

- A runtime check to ensure that all atoms have a nonzero radius. API kills
  program execution if this is the case.
- An API function to retrieve the areas/weights of the cavity finite elements.
  The values in the returned array are in Bohr^2. Addresses a feature request
  from @shoefener (Issue #13)
- The standalone executable `run_pcm` is now tested within the unit tests
  suite. The tests cover the cases where the cavity is given implicitly,
  explicitly or by substitution of radii on chosen atoms.

### Changed

- Boundary integral operators classes learnt to accept a scaling factor for the
  diagonal elements of the approximate collocation matrices. The change is
  reflected in the Green's funtion classes and in the input parsing. Addresses
  a feature request from @shoefener (Issue #16)
- `GePolCavity` learnt to print also the list of spheres used to generate the
  cavity.
- Different internal handling of conversion factors from Bohr to Angstrom.
- CMake minimum required version is 2.8.10
- `Atom`, `Solvent` and `Sphere` are now PODs. The radii and solvent lists are free
  functions.
- `PCMSOLVER_ERROR` kills program execution when an error arises but does not
  use C++ exceptions.
- `include`-s are now specified on a per-directory basis (see programmers'
  manual for a more detailed explanation)
- Default types for template paramters `DerivativeTraits`,  `IntegratorPolicy`
  and `ProfilePolicy` are now given for the Green's functions classes. This
  reduced the verbosity in instatiating these objects significantly.

### Known Issues

- The new printer in `GePolCavity` might not work properly when an explicit list
  of spheres is provided in the input.
- On Ubuntu 12.10, 32 bit the Intel compiler version 2013.1 produces a faulty
  library. It is possibly a bug in the implementation of `iso_c_binding`, see
  Issue #25

### Removed

- `SurfaceFunction` as a class is no longer available. We keep track of surface
  functions at the interface level _via_ a label-vector map.

## [Version 1.1.0] - 2016-02-07

### Added

- Green's function for diffuse interfaces in spherical symmetry

### Changed

- CMake minimum required version is 2.8.8 (2016-01-08)
- Documentation is now served [here](http://pcmsolver.readthedocs.org/)

## v1.0.4 - 2015-07-22 [YANKED]

## v1.0.3 - 2015-03-29 [YANKED]

## v1.0.2 - 2015-03-28 [YANKED]

## v1.0.1 - 2015-01-06 [YANKED]

## v1.0.0 - 2014-09-30 [YANKED]

[Unreleased]: https://github.com/PCMSolver/pcmsolver/compare/v1.1.5...HEAD
[Version 1.1.5]: https://github.com/PCMSolver/pcmsolver/compare/v1.1.4...v1.1.5
[Version 1.1.4]: https://github.com/PCMSolver/pcmsolver/compare/v1.1.3...v1.1.4
[Version 1.1.3]: https://github.com/PCMSolver/pcmsolver/compare/v1.1.2...v1.1.3
[Version 1.1.2]: https://github.com/PCMSolver/pcmsolver/compare/v1.1.1...v1.1.3
[Version 1.1.1]: https://github.com/PCMSolver/pcmsolver/compare/v1.1.0...v1.1.1
[Version 1.1.0]: https://github.com/PCMSolver/pcmsolver/releases/tag/v1.1.0
