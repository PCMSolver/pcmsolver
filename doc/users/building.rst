Building the module
===================

PCMSolver configuration and build process is managed through CMake.

Prerequisites and dependencies
------------------------------

A number of prerequisites and dependencies are to be satisfied to successfully build the module.
It will be here assumed that you want to perform a "full" build, i.e. you want to build the static libraries to be linked to your
QM program, the unit test suite and an offline copy of this documentation.

Compilers
~~~~~~~~~

+ a C++ compiler, compliant with the 1998 ISO C++ standard plus the 2003
  technical corrigendum and some additional defect reports.
+ a C compiler, compliant with the ISO C99 standard.
+ a Fortran compiler, compliant with the Fortran 2003 standard.

The list of primary test environments can be found in the ``README.md`` file.
It is entirely possible that using other compiler versions you might be able to
build the module. In order to ensure that you have a sane build, you will have
to run the unit test suite.

Libraries and toolchain programs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

+ CMake version 2.8.3 and higher;
+ Git version 1.7.1 and higher;
+ Python interpreter 2.4 and higher;
+ Boost libraries version 1.54.0 and higher;

.. note::

   Version 1.54.0 of Boost libraries is shipped with the module and resides in the `cmake/downloaded` subdirectory.
   Unless you want to use another version of Boost, you should not worry about satisfying this dependency.

+ `zlib <http://www.zlib.net/>`_ version 1.2 and higher (unit test suite only);
+ Doxygen version 1.7.6 and higher (documentation only)
+ Perl (documentation only)
+ PyYAML Python module (documentation only)
+ Sphinx (documentation only)

PCMSolver relies on the Eigen template libraries version 3.0.0 and higher.
Version 3.2.0 of Eigen libraries is shipped with the module and resides in the `external` subdirectory.

Configuration
-------------

Configuration is managed through the front-end script `setup` residing in the repository main directory.
Issuing:

.. code-block:: bash

   ./setup [options] [build path]

will create the build directory in build path and run CMake with the given options. By default, files are configured in
the `build` directory.
The `-h` or `--help` option will list the available options and their effect.
Options can be forwarded directly to CMake by using the `--cmake-options` flag and listing the `-D...` options.
Usually the following command is sufficient to get the configuration done for a debug build, including
compilation of the unit test suite:

.. code-block:: bash

   ./setup --type=debug

The unit tests suite is **always** compiled in standalone mode, unless the `-DENABLE_TESTS=OFF` option is
forwarded to CMake.

Getting Boost
~~~~~~~~~~~~~

You can get Boost libraries in two ways:

 + already packaged by your Linux distribution or through MacPorts/Brew;
 + by downloading the archive from http://www.boost.org/ and building it yourself.

In case your distribution packages a version older than 1.54.0 you might chose
to either build Boost on your own or to rely on the automated build of the
necessary Boost libraries when compiling the module (recommended).  Full
documentation on how to build Boost on Unix variants is available
`here <http://www.boost.org/doc/libs/1_56_0/more/getting_started/unix-variants.html>`_.
It is here assumed that the user **does not** have root access to the machine
and will install the libraries to a local prefix, a subdirectory of
`/home/user-name` tipically.
Once you've downloaded and unpacked the archive, run the bootstrap script to configure:

.. code-block:: bash

   cd path/to/boost
   ./bootstrap.sh --prefix=/home/user-name/boost

Running `./bootstrap.sh --help` will list the available options for the script. To build run:

.. code-block:: bash

   ./b2 install

This might take a while. After a successful build you will find the headers in
`/home/user-name/boost/include` and libraries in `/home/user-name/boost/lib`
Now, you will have Boost in a nonstandard location. Without hints CMake will
not be able to find it and configuration of `PCMSolver` will fail.  To avoid
this, you will have to pass the location of the headers and libraries to the
setup script, either with:

.. code-block:: bash

   ./setup --boost-headers=/home/user-name/boost/include --boost-libs=/home/user-name/boost/lib

or with:

.. code-block:: bash

   ./setup -DBOOST_INCLUDEDIR=/home/user-name/boost/include -DBOOST_LIBRARYDIR=/home/user-name/boost/lib

Advanced configuration options
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

As of version 3.1.0 Eigen can fall back to certain BLAS/LAPACK routines as
implemented in Intel's Math Kernel Library (MKL).  You can enable this feature
passing the option `-DENABLE_EIGEN_MKL=ON` to the setup script. You need at least
version 10.3 of Intel MKL.  Be sure to have read
the relevant `guide <http://eigen.tuxfamily.org/dox/TopicUsingIntelMKL.html>`_

.. warning::

   By itself, the module **does not** need to link to any mathematical
   library. This option makes sense only when linking the module to a QM program
   that links to Intel MKL. **Use at your own risk!** There is no proof that
   enabling this option leads to significant speedups!

Build and test
--------------

To compile and link, just go to the build directory and run:

.. code-block:: bash

   make -j N

where `N` is the number of cores you want to use when building.

.. note::

   Building on more than one core can sometimes result in a "race condition"
   and a crash. If that happens, please report the problem as an issue on our
   issue tracker on GitHub. Running `make` on a single core might get you through
   compilation.

To run the whole test suite:

.. code-block:: bash

   ctest -j N

You can also use CTest to run a specific test or a set of tests. For example:

.. code-block:: bash

   ctest -R gepol

will run all the test containing the string "gepol" in their name.

If Doxygen was found, an offline copy of this documentation can be built by:

.. code-block:: bash

   make doc

and visualized by opening the `doc/html/index.html` file in your browser.
