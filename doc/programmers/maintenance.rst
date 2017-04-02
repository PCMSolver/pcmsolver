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

Updating the copyright notice
-----------------------------

You need to have access to the ``pcmsolvermeta`` repository to update the
copyright notice.  The copyright notice text is in the file
``copyright_notice.txt`` in the ``tools`` directory.  The script
``update_copyright.py`` will extract the text from the file, create the
appropriate header and perform the update on the files in the subdirectory
where it is invoked.

.. warning::
   The copyright notice on top of the Config.hpp.in file needs to be **manually** updated!

Release process
---------------

We have two repositories one public for the release, hosted on `GitHub
<https://github.com/PCMSolver/pcmsolver>`_ and one private for the
development, hosted on `GitLab <https://gitlab.com/PCMSolver/pcmsolver>`_.
At release time the master branch on the private repository is synced to that
of the public repository.

.. warning::
   This means that **WHATEVER** is on master at release time is considered
   ready for release.  Protection of functionality happens **EXCLUSIVELY** by
   making use of branches/forks on the private repository.

You need to compile the to-be-released code and run the unit test suite.  If
compilation works and all unit tests are passing then the code is ready to be
released:

.. code-block:: bash

   git push Origin release

Notice that ``Origin`` has been spelled with a capital ``O`` the reason being
that the release branch gets pushed both to the private and the public
repositories (`trick explanation
<http://stackoverflow.com/questions/849308/pull-push-from-multiple-remote-locations>`_)
In brief, you need to have a ``.git/config`` file that resembles the following:

.. code-block:: bash

   [remote "origin"]
       url = git@gitlab.com:PCMSolver/pcmsolver.git
       fetch = +refs/heads/*:refs/remotes/origin/*
   [remote "GitHub"]
       url = git@github.com:PCMSolver/pcmsolver.git
       fetch = +refs/heads/*:refs/remotes/GitHub/*
   [remote "Origin"]
       url = git@gitlab.com:PCMSolver/pcmsolver.git
       url = git@github.com:PCMSolver/pcmsolver.git
