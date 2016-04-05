# Change Log

## [Unreleased]

## [v1.1.1] (2016-03-10)

### Added

- A runtime check to ensure that all atoms have a nonzero radius.
API kills program execution if this is the case.
- An API function to retrieve the areas/weights of the cavity finite elements.
The values in the returned array are in Bohr^2. Addresses a feature request from @shofener (Issue #13)
- The standalone executable `run_pcm` is now tested within the unit tests suite.
The tests cover the cases where the cavity is given implicitly, explicitly or by substitution of
radii on chosen atoms.

### Changed

- Boundary integral operators classes learnt to accept a scaling factor for the diagonal
elements of the approximate collocation matrices. The change is reflected in the
Green's funtion classes and in the input parsing. Addresses a feature request from @shofener (Issue #16)
- GePolCavity learnt to print also the list of spheres used to generate
the cavity.
- Different internal handling of conversion factors from Bohr to Angstrom.
- CMake minimum required version is 2.8.10
- Atom, Solvent and Sphere are now PODs. The radii and solvent lists are free functions.
- `PCMSOLVER_ERROR` kills program execution when an error arises but does
not use C++ exceptions.
- `include`-s are now specified on a per-directory basis (see programmers' manual
for a more detailed explanation)
- Default types for template paramters `DerivativeTraits`,  `IntegratorPolicy` and `ProfilePolicy`
are now given for the Green's functions classes. This reduced the verbosity in instatiating
these objects significantly.

### Known Issues

- The new printer in GePolCavity might not work properly when an explicit list
of spheres is provided in the input.
- On Ubuntu 12.10, 32 bit the Intel compiler version 2013.1 produces a faulty library.
It is possibly a bug in the implementation of `iso_c_binding`, see Issue #25

### Removed

- SurfaceFunction as a class is no longer available. We keep track of surface functions
at the interface level _via_ a label-vector map.

## [v1.1.0] (2016-02-07)

### Added

- Green's function for diffuse interfaces in spherical symmetry

### Changed

- CMake minimum required version is 2.8.8 (2016-01-08)
- Documentation is now served [here](http://pcmsolver.readthedocs.org/)

## v1.0.4 (2015-07-22) [YANKED]

## v1.0.3 (2015-03-29) [YANKED]

## v1.0.2 (2015-03-28) [YANKED]

## v1.0.1 (2015-01-06) [YANKED]

## v1.0.0 (2014-09-30) [YANKED]

[Unreleased]: https://github.com/PCMSolver/pcmsolver/compare/v1.1.1...HEAD
[v1.1.1]: https://github.com/PCMSolver/pcmsolver/compare/v1.1.0...v1.1.1
[v1.1.0]: https://github.com/PCMSolver/pcmsolver/releases/tag/v1.1.0
