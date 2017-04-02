Maintenance
===========

Description and how-to for maintenance operations.
Some of the maintenance scripts have been moved to the `pcmsolvermeta
repository <https://gitlab.com/PCMSolver/pcmsolvermeta>`_

Bump version
------------

Version numbering follows the guidelines of `semantic versioning <http://semver.org/>`_
To update, change the relevant field in the ``README.md`` file.

Changelog
---------

We follow the guidelines of `Keep a CHANGELOG <http://keepachangelog.com/>`_
On all **but** the release branches, there is an ``Unreleased`` section
under which new additions should be listed.
To simplify perusal of the ``CHANGELOG.md``, use the following subsections:

1. ``Added`` for new features.
2. ``Changed`` for changes in existing functionality.
3. ``Deprecated`` for once-stable features removed in upcoming releases.
4. ``Removed`` for deprecated features removed in this release.
5. ``Fixed`` for any bug fixes.
6. ``Security`` to invite users to upgrade in case of vulnerabilities.

Updating Eigen distribution
---------------------------

The C++ linear algebra library Eigen comes bundled with the module. To update
the distributed version one has to:

1. download the desired version of the library to a scratch location. Eigen's
   website is: http://eigen.tuxfamily.org/
2. unpack the downloaded archive;
3. go into the newly created directory and create a build directory;
4. go into the newly created build directory and type the following (remember
   to substitute @PROJECT_SOURCE_DIR@ with the actual path)

   .. code-block:: bash

	  cmake .. -DCMAKE_INSTALL_PREFIX=@PROJECT_SOURCE_DIR@/external/eigen3

Remember to commit and push your modifications.

Release process
---------------

.. warning::
   **Incomplete or outdated information!**

Releases in a `X.Y.Z` series are annotated tags on the corresponding branch.
