# Change Log

## [Unreleased]

## [Version 1.2.3] - 2019-02-14

### Added

- A `pcmsolver_citation` function was added which accepts a function pointer to
  printing facilities used in the host program. `pcmsolver_citation` will print
  versioning and citation information for the library. Please, include a call to
  this function and print the returned string somewhere in your code, so that
  users are aware of what literature to cite when using the library.

### Changed

- The list of finite elements return from the PEDRA Fortran code is now pruned
  to remove finite elements with area less than 1.0e-4 This avoids numerical
  artifacts in the formation of the PCM matrices.
  The number of pruned finite elements is reported when printing the cavity.
- The `pcmsolver_print` function now only prints out the set up for the PCM calculation.

## [Version 1.2.2] - 2019-01-27

### Changed

- Perform Gauss' theorem check on nuclear ASC upon initialization. In case a
  discrepancy larger than 10^-3 is detected in the absolute difference between the
  Gauss' theorem and computed values, calculation is aborted as this points to
  numerical problems, most likely in the construction of the cavity.
- Add check and test for S matrix positive-definiteness. Certain cavity
  discretizations lead to the S matrix not being symmetric positive-definite. We
  abort if that's the case.

## [Version 1.2.1] - 2018-05-02

### Fixed

- Import of local copy of the `pyparsing` module in Getkw. Thanks @loriab for
  the fix.
- The `conf.py` documentation build script is restored to normal operation.

## [Version 1.2.0] - 2018-04-27

### Added

- Green’s function for a spherical nanoparticle with **real** permittivity.
  The Green's function is known in analytical form from the work of
  Messina (J. Chem. Phys. 2002, 117 (24), 11062) and
  Delgado _et al._ (J. Chem. Phys. 2013, 139 (2), 024105)

### Deprecated

- C++03 support is effectively **deprecated** in favor of C++11. Support for GCC
  4.6 and earlier has been dropped. Please consider upgrading your C++ compiler
  to a [fully standard-compliant one](https://en.cppreference.com/w/cpp/compiler_support#cpp11)

### Changed

- The versioning machinery has been updated. The update was inspired by the
  `versioner.py` tool devised by @loriab.
- The documentation building scripts `conf.py` and `cloc_tools.py` have been
  thoroughly refactored and simplified:
  * `cloc_tools.py` is now called `cloc_wrapper.py`. It is no longer configured
  and can be found in the `cloc_tools` subfolder of `doc`.
  * A local build of documentation will only work if run from within the `doc` folder:
  ```
  sphinx-build . _build
  ```
  This choice was made to simplify the set up of the ReadTheDocs and local
  documentation building procedures and to minimize the chances of breaking
  either.
- The Fortran API bindings file `pcmsolver.f90` is now installed alongside the
  `pcmsolver.h` header file. The users will have to compile it explicitly to
  get the type checking from the API redeclaration in Fortran 90.
  The file is always installed.
- The `ENABLE_Fortran_API` configuration option has been renamed
  `TEST_Fortran_API`, since it now only triggers compilation of the
  `Fortran_host` test case.
- **BREAKING CHANGE** The minimum required version of CMake is now 3.3

### Fixed

- `std::string`-s are now used in the `Meddle` object functions manipulating
  the surface functions map. Explicit casts from `const char *` to `std::string`
  are handled in the API functions.
- Properly enforce `const`-correctness of the `Meddle` object and of its usage
  in the context API.
- [Cholesky decomposition](http://eigen.tuxfamily.org/dox/classEigen_1_1LLT.html) is used
  in CPCMSolver to get the inverse of the S matrix. The robust Cholesky (LDLT)
  previously used is broken with the latest version of the Intel compilers.

## [Version 1.2.0-rc1] - 2018-03-02

### Added

- Double logarithmic scale for the integration of spherical diffuse
  interfaces: much more stable than the previous version, allowing for
  Runge-Kutta 4 integrator.
- A new CMake module `options_wrappers.cmake` that adds new wrapper macros for
  the CMake `option` command.

### Fixed

- Bug in the diffuse interface Green's function. Contrary to the sharp
  interface case, it is wrong to remove the monopole, which becomes
  identically zero when the corresponding differential equation is
  solved in extreme cases (e.g. charge far away from the sphere).
- Visibility of symbols in the shared library is _finally_ handled properly.
  The necessary flags to the C++ compiler were not set for the subtargets built
  as CMake `OBJECT` libraries. This results in a modest decrease in library
  size.
- `RPATH` handling for the standalone executable `run_pcm`.
  For the executable in `<build_dir>/bin` the `RPATH` will contain
  `<build_dir>/lib64` (or `<build_dir>/lib`) as path to `libpcm.so.1`.
  For the executable in `<install_prefix>/bin` the `RPATH` will contain
  `<install_prefix>/lib64` (or `<install_prefix>/lib`) as path to `libpcm.so.1`.
- Code coverage analysis was restored. We now use
  [Codecov](https://codecov.io). Thanks @arnfinn :tada:
- **BREAKING CHANGE** The layout of the installation for the Python scripts has
  been changed. This addresses issue #116. Given the user-defined installation
  prefix (`<prefix>`) one would obtain:
  ```
  <prefix>/pcmsolver
           ├── bin
           │   ├── go_pcm.py
           │   ├── plot_cavity.py
           │   └── run_pcm
           ├── include
           │   └── PCMSolver
           │       ├── bi_operators
           │       ├── cavity
           │       ├── Citation.hpp
           │       ├── Config.hpp
           │       ├── Cxx11Workarounds.hpp
           │       ├── ErrorHandling.hpp
           │       ├── external
           │       ├── GitInfo.hpp
           │       ├── green
           │       ├── interface
           │       ├── LoggerInterface.hpp
           │       ├── PCMInput.h
           │       ├── PCMSolverExport.h
           │       ├── pcmsolver.h
           │       ├── PhysicalConstants.hpp
           │       ├── solver
           │       ├── STLUtils.hpp
           │       ├── TimerInterface.hpp
           │       └── utils
           ├── lib64
           │   ├── libpcm.a
           │   ├── libpcm.so -> libpcm.so.1
           │   ├── libpcm.so.1
           │   └── python
           │       └── pcmsolver
           └── share
               └── cmake
                   └── PCMSolver
  ```

### Changed

- As a result of the visibility change, unit tests can only be linked against
  the static library, since all symbols are always visible in a static archive
  library.
- **BREAKING CHANGE** The `pcmsolver.py` script/module **was removed** to
  address issues #111 and #112. The `main` portion of the script has now been
  separated into a `go_pcm.py` script which can be used in the same way as
  `pcmsolver.py` was used before: to parse the PCMSolver input and to run the
  module standalone. This separation greatly simplifies the use as a module
  within the Python launcher scripts of host programs.

## [Version 1.1.12] - 2018-01-20

### Added

- The function `pcmsolver_fstring_to_carray` was added to the Fortran bindings
  for the library. As the name suggests, this function "translates" a Fortran
  string into a C `char` array. This supersedes and replaces the
  `pcmsolver_f2c_string` function.
- An updated version of the `pcmsolver_new` API function has been introduced.
  The full name of the function is `pcmsolver_new_v1112` (to avoid API breakage)
  and accepts as additional argument the name of the parsed input file
  generated by `pcmsolver.py`. This allows to decoupled parsing and writing out
  the parsed file, in case directory-specific operations have to be performed
  by the host program.

### Changed

- Dispatching of Green’s functions in the factory is purely label-based. The
  input reader generates a string that holds the type of the Green’s function,
  the strategy for calculating its normal derivatives and the dielectric
  profile (if needed): `TYPE_DERIVATIVE_PROFILE`. The Green’s functions are
  subscribed to the factory using this label, which is then also the one used
  for their creation. Iteration over _lists of types_ to dispatch the correct
  template parameters is thus completely avoided. This removes the dependency
  on the metafunctions in `ForId.hpp` and the dependency on Boost.MPL, the
  Boost metaprogramming library.
- Dispatching of quadrature rules for the numerical boundary integral operator
  has been rewritten, removing the dependency on Boost.MPL.
- The input wrapping data `struct`-s have been modified to include the string
  defining the type of object and to remove now unused information.
- Documentation building is fully handled _via_ `sphinx-build`: CMake will no longer generate a `doc` build target.
- Simplified `.travis.yml` and got rid of Conda to handle multiple Python versions.
- Moved python scripts from `bin/` to `bin/pcmsolver/` and introduced `bin/pcmsolver/__init__.py`.
- The `parse_pcm_input` function returns the parsed input as a string without
  saving it to file. The old behavior can be recovered by setting the
  `write_out` argument to `True`.
- The [Catch unit test framework](https://github.com/philsquared/Catch) has
  been updated to its latest version still supporting C++03
  [v1.11.0](https://github.com/philsquared/Catch/releases/tag/v1.11.0)
- Travis continuous integration is _no longer_ run on Mac infrastructure.

### Fixed

- Documentation building on ReadTheDocs is fully functional again, thanks @arnfinn :tada:
  The build had been failing for a while since docs were generated for all files, including
  documentation files from previous build. Besides, source code doxygen blocks were not
  extracted when inside namespaces.
- Compiler warning in the `Sphere` class due to a redeclaration of `operator<<`.
- The `plot_cavity.py` script is now Python 3 compatible.

### Removed

- `FortranCUtils` header and source files were removed. The `pcmsolver_f2c_string` and `pcmsolver_c2f_string`
  functions are thus gone. For the former, use the replacement
  `pcmsolver_fstring_to_carray` function provided in the Fortran bindings to the library.
- Some unused files have been removed:
  * `Interpolation.hpp`
  * `Interpolation.cpp`
  * `Vector2.hpp`
  * `Vector3.hpp`
  * `ForId.hpp`
  * `DerivativeUtils.hpp`
- The `cavityType()`, `solverType()`, `greenInsideType()`, `greenOutsideType()`
  and `integratorType()` functions in the `Input` object have been removed.
  This information is now wrapped into the corresponding input-wrapping
  `struct`-s: `CavityData`, `SolverData`, `GreenData` and `BIOperatorData`.

### Deprecated

- We are in the process of removing the dependency on the Boost libraries. The
  only dependency on Boost allowed in the 1.Y.Z release series will be for
  replacing C++11 functionality that is missing from the C++98 standard. For
  such code, we will require a pure C++11 counterpart.

## [Version 1.1.11] - 2017-10-25

### Added

- A Python script, using Matplotlib, to plot the cavity.
  The script can also color-map the finite elements according to the values of
  a surface function.
- The input learnt to parse the additional `ChargeDistribution` section.
  It is possible to specify a classical charge distribution of point multipoles.
  This can be an additional source of electrostatic potential for the calculation
  of the ASC.
- Restored compilation for g++ < v5.1.
- [Ninja](https://ninja-build.org/) can be used as a generator.
  Notice that at least [CMake 3.7.2](https://cmake.org/cmake/help/v3.7/generator/Ninja.html#fortran-support)
  **and** the [Kitware-maintained](https://github.com/Kitware/ninja) version of
  Ninja are required to successfully compile.

### Changed

- Use [`#pragma once`](https://en.wikipedia.org/wiki/Pragma_once) instead of
  `#ifndef, #define, #endif` to guard against multiple inclusion of header files.
- The uppercased contents of the `.pcm` input file are written to a temporary
  file, instead of overwriting the user provided file. The temporary file is
  removed after it has been parsed. Fixes #91 as noted by @ilfreddy.
- Use Runge-Kutta-Fehlberg 7(8) ODE solver to integrate the radial equation
  in the spherical diffuse Green's function class.

### Fixed

- A bug in the initialization of a restart cavity from old `.npz` files.
  Currently, the `.npz` file saves sphere center, arcs and vertices of each
  finite element. This information is in fact needed to plot the cavity using
  the Python script in `tools`. Older `.npz` files did not contain this
  information and could not be read in. The additional information is read in
  as arrays of zeros in case it is not present on file.

## [Version 1.1.10] - 2017-03-27

### Changed

- Updated the `cloc.pl` script to version 1.72
- Simplified the internal structure of the `Meddle` and `Input` objects.
- Export dependency on Zlib for the static libraries. Thanks @loriab for the pull request
  fixing [a build problem within Psi4](http://forum.psicode.org/t/crc32-undefined-symbol-at-runtime-when-built-with-pcmsolver-gcc-4-9-4/449/7)

## [Version 1.1.9] - 2017-02-16

### Changed

- PCMSolver is now exported as a proper [CMake target](https://cmake.org/cmake/help/v3.0/manual/cmake-buildsystem.7.html)
  See PR #38 for details. Thanks @loriab for the work.
- The Python scripts shipped with the library are now Python 2 and Python 3 compatible.
- `Factory` is no longer implemented as a Singleton.
- The [Catch unit test framework](https://github.com/philsquared/Catch) has
  been updated to its latest version
  [v1.7.2](https://github.com/philsquared/Catch/releases/tag/v1.7.2)
- Updated the version of Eigen bundled with the code.
  The minimum required version of Eigen is still 3.3.0, but we ship
  [Eigen 3.3.2](http://eigen.tuxfamily.org/index.php?title=ChangeLog#Eigen_3.3.2)

### Fixed

- Revert to use [Robust Cholesky decomposition](https://eigen.tuxfamily.org/dox/classEigen_1_1LDLT.html)
  to compute the inverse of the S matrix in `CPCMSolver`.

## [Version 1.1.8] - 2017-02-06

### Added

- Namespaces for all of the internal code have been introduced.
  The top-level namespace is `pcm`. At finer levels the namespaces have the same
  names as the respective subdirectories. Read the programmers' documentation
  for further details on the use of namespaces in the project.
  Related to issue #34 on [GitHub] and #60 on [GitLab].
- The top-level, convenience header `BoundaryIntegralOperator.hpp` includes all
  subclasses and utility headers in the `bi_operators` subdirectory.
  Related to issue #34 on [GitHub] and #60 on [GitLab].
- The top-level, convenience header `Cavity.hpp` includes all subclasses and
  utility headers in the `cavity` subdirectory.
  Related to issue #34 on [GitHub] and #60 on [GitLab].
- The top-level, convenience header `Green.hpp` includes all subclasses and
  utility headers in the `green` subdirectory.
  Related to issue #34 on [GitHub] and #60 on [GitLab].
- The top-level, convenience header `Solver.hpp` includes all subclasses and
  utility headers in the `solver` subdirectory.
  Related to issue #34 on [GitHub] and #60 on [GitLab].

### Changed

- The abstract base class for the boundary integral operator integrators has
  been renamed `IBoundaryIntegralOperator`.
  The relevant factory is bootstrapped upon creation of the `Meddle` object,
  i.e. at library initialization.
  Related to issue #34 on [GitHub] and #60 on [GitLab].
- The abstract base class for the cavities has been renamed `ICavity`.
  The relevant factory is bootstrapped upon creation of the `Meddle` object,
  i.e. at library initialization.
  Related to issue #34 on [GitHub] and #60 on [GitLab].
- The abstract base class for the solvers has been renamed `ISolver`.
  The relevant factory is bootstrapped upon creation of the `Meddle` object,
  i.e. at library initialization.
  Related to issue #34 on [GitHub] and #60 on [GitLab].
- The `typedef` for numerical differentiation in the Green's function classes
  has been renamed `Stencil` to avoid name clashes with the `Numerical`
  boundary integral operator type.

### Fixed

- A bug in the selection of the extended diagnostics flags for the GNU C++
  compiler. These flags are now enabled only for versions >= 5.1.0 and when the
  C++11 standard is enable. Fixes issue #36 on [GitHub] and #62 on [GitLab].
- A bug in the initialization of the factory for the cavity classes was fixed.
  The bug manifested only in the static library `libpcm.a`
  Fixes issue #34 on [GitHub] and #60 on [GitLab].

## [Version 1.1.7] - 2016-12-01

### Added

- A pre-commit hook in `.githooks/pre-commit-clang-format` checking that the
  format of C++ header and source files conforms to the project style. The hook
  uses `clang-format` and the style mandated by the `.clang-format` file to
  check files in the tree. Commit is rejected if one or more files are non
  compliant. The hook generates a patch and shows the command needed to apply
  it.
  To enable the hooks you need to have a `.git/hooks/pre-commit` file
  containing this line `.githooks/pre-commit`
  _NOT recommended_ The hook can be skipped by passing the `--no-verify` option to `git commit`
- A pre-commit hook in `.githooks/pre-commit-license-maintainer` checking the
  license headers. The hook is configured based on the `.gitattributes` file.
  The hook will check the license headers and amend them, either by updating
  the year and authors information or by adding the whole header.
  To enable the hooks you need to have a `.git/hooks/pre-commit` file
  containing this line `.githooks/pre-commit`
  _NOT recommended_ The hook can be skipped by passing the `--no-verify` option to `git commit`
- An `UNUSED` preprocessor macro to mark arguments as unused.
- An `UNUSED_FUNCTION` preprocessor macro to mark functions as unused.
- A set of preprocessor macros with Git information (`GIT_COMMIT_HASH`,
  `GIT_COMMIT_AUTHOR`, `GIT_COMMIT_DATE` and `GIT_BRANCH`) is automatically
  generated and saved in the `git_info.h` header file.
- An API function to print the contents of a surface function to the host
  program output.
- An API function to get the dipole moment, relative to the origin, due to the ASC
  on the cavity. Both the norm and the components can be obtained.
- `.gitattributes` now instructs Git to ignore binary files in diff operations.
  PNG files are diff-ed using EXIF information. To set this up properly,
  install an EXIF tool on your machine and run `git config diff.exif.textconv
  exiftool` in your local copy of the repository.

### Changed

- The Fortran bindings file has been renamed `pcmsolver.f90`.
- The Green's function, solver and boundary integral operator classes have been
  radically redesigned. This avoids coupling between integrators and Green's
  function that existed in the previous design.
  See the [Green's function code
  reference](http://pcmsolver.readthedocs.io/en/latest/code-reference/greens-functions.html)
  for a more detailed explanation.
- **BREAKING CHANGE** The minimum required version of Eigen is now 3.3.0
  The version bundled with the code has been accordingly updated.
- The `PCMSOLVER_ERROR` macro now takes only one argument and prints out a more
  informative error message.
- Switched to the latest version of
  [Autocmake](http://autocmake.readthedocs.io/) The configuration file is now YAML-based.
  The PyYAML module is thus required.
- The extended diagnostic flags `-Wsuggest-attribute=pure
  -Wsuggest-attribute=const -Wsuggest-attribute=noreturn -Wsuggest-final-types
  -Wsuggest-final-methods -Wsuggest-override -Wuseless-cast
  -Wunsafe-loop-optimizations` are always set when using the GNU C++ compiler
  in a debug configuration.
- The C++11 compatibility CMake macros now check for the availability of the
  `noreturn` attribute. A workaround macro, accessible _via_ `__noreturn`, has
  been added to the `Cxx11Workarounds.hpp` header file.
- **BREAKING CHANGE** The ouput flushing function must be passed explicitly as
  a function pointer to the `pcmsolver_new` function during library
  initialization.
  The function pointer has the signature
  `typedef void (*HostWriter)(const char * message)`
  thus accepting a single argument instead of the previous two.
- [GNU standard installation
  directories](http://www.gnu.org/prep/standards/html_node/Directory-Variables.html)
  have been imposed, thanks to work by @loriab.
  Given a prefix, header files are now installed to `include/pcmsolver`,
  executables to `bin`, libraries to `lib` and scripting tools to `share`.
  The install prefix and the installation directories can be specified by the
  `--prefix`, `--bindir`, `--libdir`, `--includedir` and `--datadir` options to
  the `setup.py` script (or the corresponding CMake variables)

## [Version 1.1.6] - 2016-09-20

### Added

- A function returning a molecule object for the water molecule.

### Changed

- [Cholesky decomposition](http://eigen.tuxfamily.org/dox/classEigen_1_1LLT.html) is used
  whenever the inverse of the S matrix has to be calculated.
  The S matrix is self-adjoint, positive-definite and the LLT decomposition is
  faster than LDLT.

### Fixed

- Some inconsistencies in input reading from host and a related memory leak in the radii
  initialization.

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

[Unreleased]: https://github.com/PCMSolver/pcmsolver/compare/v1.2.3...HEAD
[Version 1.2.3]: https://github.com/PCMSolver/pcmsolver/compare/v1.2.2...v1.2.3
[Version 1.2.2]: https://github.com/PCMSolver/pcmsolver/compare/v1.2.1...v1.2.2
[Version 1.2.1]: https://github.com/PCMSolver/pcmsolver/compare/v1.2.0...v1.2.1
[Version 1.2.0]: https://github.com/PCMSolver/pcmsolver/compare/v1.2.0-rc1...v1.2.0
[Version 1.2.0-rc1]: https://github.com/PCMSolver/pcmsolver/compare/v1.1.12...v1.2.0-rc1
[Version 1.1.12]: https://github.com/PCMSolver/pcmsolver/compare/v1.1.11...v1.1.12
[Version 1.1.11]: https://github.com/PCMSolver/pcmsolver/compare/v1.1.10...v1.1.11
[Version 1.1.10]: https://github.com/PCMSolver/pcmsolver/compare/v1.1.9...v1.1.10
[Version 1.1.9]: https://github.com/PCMSolver/pcmsolver/compare/v1.1.8...v1.1.9
[Version 1.1.8]: https://github.com/PCMSolver/pcmsolver/compare/v1.1.7...v1.1.8
[Version 1.1.7]: https://github.com/PCMSolver/pcmsolver/compare/v1.1.6...v1.1.7
[Version 1.1.6]: https://github.com/PCMSolver/pcmsolver/compare/v1.1.5...v1.1.6
[Version 1.1.5]: https://github.com/PCMSolver/pcmsolver/compare/v1.1.4...v1.1.5
[Version 1.1.4]: https://github.com/PCMSolver/pcmsolver/compare/v1.1.3...v1.1.4
[Version 1.1.3]: https://github.com/PCMSolver/pcmsolver/compare/v1.1.2...v1.1.3
[Version 1.1.2]: https://github.com/PCMSolver/pcmsolver/compare/v1.1.1...v1.1.3
[Version 1.1.1]: https://github.com/PCMSolver/pcmsolver/compare/v1.1.0...v1.1.1
[Version 1.1.0]: https://github.com/PCMSolver/pcmsolver/releases/tag/v1.1.0

[GitHub]: https://github.com/PCMSolver/pcmsolver
[GitLab]: https://gitlab.com/PCMSolver/pcmsolver
