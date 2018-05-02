Versioning and minting a new release
====================================

Our versioning machinery is based on a modified version of the ``versioner.py``
script devised by Lori A. Burns (Georgia Tech) for the `Psi4
<http://www.psicode.org>`_ quantum chemistry code.
The documentation that follows is also adapted from the corresponding Psi4
documentation, available at `this link <http://www.psicode.org/psi4manual/1.1/manage_git.html>`_

This guide will walk you through the actions to perform to mint a new release
of the code. Version numbering follows the guidelines of `semantic versioning
<http://semver.org/>`_. The allowed format is ``MAJOR.MINOR.PATCH-DESCRIBE``,
where ``DESCRIBE`` can be a string describing a prerelease state, such as
``rc2``, ``alpha1``, ``beta3`` and so forth.

Minting a new release
---------------------

The ``tools/metadata.py`` file records the versioning information for the current
release. The information in this file is used by the ``versioner.py`` script to
compute a *unique version number* for development snapshots.

.. note::

   To correctly mint a new release, you will have to be on the latest release
   branch of (i) a direct clone or (ii) clone-of-fork with release branch
   up-to-date with upstream (including tags!!!) and with upstream as remote.

This is the step-by-step guide to releasing a new version of PCMSolver:

#. **DECIDE** an upcoming version number, say ``1.2.0``.
#. **TIDY UP** ``CHANGELOG.md``:

   * **SET** the topmost header to the upcoming version number and release date.

     ::

       ## [Version 1.2.0] - 2018-03-31

   * **CHECK** that the links at the bottom of the document are correct.

     ::

       [Unreleased]: https://github.com/PCMSolver/pcmsolver/compare/v1.2.0...HEAD
       [Version 1.2.0]: https://github.com/PCMSolver/pcmsolver/compare/v1.2.0-rc1...v1.2.0
       [Version 1.2.0-rc1]: https://github.com/PCMSolver/pcmsolver/compare/v1.1.12...v1.2.0-rc1

#. **UPDATE** the ``AUTHORS.md`` file:

   * Run ``git shortlog -sn`` and cross-check with the current contents of ``AUTHORS.md``.
     Edit where necessary and don't forget to include, where
     available, the GitHub handle. Authors are ordered by the number of commits.
   * Update the revision date at the bottom of this file.

     ::

       >>> cat AUTHORS.md
       ## Individual Contributors

       - Roberto Di Remigio (@robertodr)
       - Luca Frediani (@ilfreddy)
       - Monica Bugeanu (@mbugeanu)
       - Arnfinn Hykkerud Steindal (@arnfinn)
       - Radovan Bast (@bast)
       - T. Daniel Crawford (@lothian)
       - Krzysztof Mozgawa
       - Lori A. Burns (@loriab)
       - Ville Weijo (@vweijo)
       - Ward Poelmans (@wpoely86)

       This list was obtained 2018-03-02 by running `git shortlog -sn`

#. **CHECK** that the ``.mailmap`` file is up-to-date.
#. **CHECK** that the documentation builds locally.
#. **ACT** to check all the changed files in.
#. **OBSERVE** current versioning state

   * https://github.com/PCMSolver/pcmsolver/releases says ``v1.2.0-rc1`` & ``9a8c391``

    ::

      >>> git tag
      v1.1.0
      v1.1.1
      v1.1.10
      v1.1.11
      v1.1.12
      v1.1.2
      v1.1.3
      v1.1.4
      v1.1.5
      v1.1.6
      v1.1.7
      v1.1.8
      v1.1.9
      v1.2.0-rc1

      >>> cat tools/metadata.py
      __version__ = '1.2.0-rc1'
      __version_long = '1.2.0-rc1+9a8c391'
      __version_upcoming_annotated_v_tag = '1.2.0'
      __version_most_recent_release = '1.1.12'


      def version_formatter(dummy):
          return '(inplace)'

      >>> git describe --abbrev=7 --long --always HEAD
      v1.2.0-rc1-14-gfc02d9d

      >>> git describe --abbrev=7 --long --dirty
      v1.2.0-rc1-14-gfc02d9d-dirty

      >>> python tools/versioner.py
      Defining development snapshot version: 1.2.0.dev14+fc02d9d (computed)
      1.2.0.dev14 {versioning-script} fc02d9d 1.1.12.999 dirty  1.1.12 <-- 1.2.0.dev14+fc02d9d

      >>> git diff

   * Observe that current latest tag matches metadata script and git
     describe, that GH releases matches metadata script, that upcoming in
     metadata script matches current ``versioner.py`` version.

#. **ACT** to bump tag in code. The current tag is ``v1.2.0-rc1``, the imminent tag is ``v1.2.0``.

   * Edit current & prospective tag in ``tools/metadata.py``. Use your
     decided-upon tag ``v1.2.0`` and a speculative next tag, say ``v1.3.0``,
     and use 7 "z"s for the part you can't predict.

     ::

       >>> vim tools/metadata.py

       >>> git diff
       diff --git a/tools/metadata.py b/tools/metadata.py
       index 5d87b55..6cbc05e 100644
       --- a/tools/metadata.py
       +++ b/tools/metadata.py
       @@ -1,6 +1,6 @@
       -__version__ = '1.2.0-rc1'
       -__version_long = '1.2.0-rc1+9a8c391'
       -__version_upcoming_annotated_v_tag = '1.2.0'
       -__version_most_recent_release = '1.1.12'
       +__version__ = '1.2.0'
       +__version_long = '1.2.0+zzzzzzz'
       +__version_upcoming_annotated_v_tag = '1.3.0'
       +__version_most_recent_release = '1.2.0'

   * **COMMIT** changes to ``tools/metadata.py``.

     ::

       >>> git add tools/metadata.py
       >>> git commit -m "Bump version to v1.2.0"

#. **OBSERVE** undefined version state. Note the 7-character git hash for the new commit, here ``fc02d9d``.

   ::

     >>> git describe --abbrev=7 --long --always HEAD
     v1.2.0-rc1-14-gfc02d9d

     >>> git describe --abbrev=7 --long --dirty
     v1.2.0-rc1-14-gfc02d9d-dirty

     >>> python tools/versioner.py
     Undefining version for irreconcilable tags: 1.2.0-rc1 (computed) vs 1.2.0 (recorded)
     undefined {versioning-script} fc02d9d 1.2.0.999 dirty  1.2 <-- undefined+fc02d9d

#. **ACT** to bump tag in git, then bump git tag in code.

   * Use the decided-upon tag ``v1.2.0`` and the observed hash ``fc02d9d`` to
     mint a new *annotated* tag, minding that "v"s are present here.

   * Use the observed hash to edit ``tools/metadata.py`` and commit immediately.

   ::

     >>> git tag -a v1.2.0 fc02d9d -m "Version 1.2.0 released"

     >>> vim tools/metadata.py

     >>> git diff
     diff --git a/tools/metadata.py b/tools/metadata.py
     index 6cbc05e..fdc202e 100644
     --- a/tools/metadata.py
     +++ b/tools/metadata.py
     @@ -1,5 +1,5 @@
      __version__ = '1.2.0'
     -__version_long = '1.2.0+zzzzzzz'
     +__version_long = '1.2.0+fc02d9d'
      __version_upcoming_annotated_v_tag = '1.3.0'
      __version_most_recent_release = '1.2.0'

     >>> python tools/versioner.py
     Amazing, this can't actually happen that git hash stored at git commit.

     >>> git add tools/metadata.py

     >>> git commit -m "Records tag for v1.2.0"

#. **OBSERVE** current versioning state. There is nothing to take note of. This
   is just a snapshot to ensure that you did not mess up.

    ::

      >>> python tools/versioner.py
      Defining development snapshot version: 1.2.0.dev1+4e0596e (computed)
      1.2.0.dev1 {master} 4e0596e 1.2.0.999   1.2 <-- 1.2.0.dev1+4e0596e

      >>> git describe --abbrev=7 --long --always HEAD
      v1.2.0-1-g4e0596e

      >>> git describe --abbrev=7 --long --dirty
      v1.2.0-1-g4e0596e

      >>> git tag
      v1.1.0
      v1.1.1
      v1.1.10
      v1.1.11
      v1.1.12
      v1.1.2
      v1.1.3
      v1.1.4
      v1.1.5
      v1.1.6
      v1.1.7
      v1.1.8
      v1.1.9
      v1.2.0-rc1
      v1.2.0

      >>> cat tools/metadata.py
      __version__ = '1.2.0'
      __version_long = '1.2.0+fc02d9d'
      __version_upcoming_annotated_v_tag = '1.3.0'
      __version_most_recent_release = '1.2.0'

      >>> cat metadata.out.py | head -8
      __version__ = '1.2.0.dev1'
      __version_branch_name = 'master'
      __version_cmake = '1.2.0.999'
      __version_is_clean = 'True'
      __version_last_release = '1.2.0'
      __version_long = '1.2.0.dev1+4e0596e'
      __version_prerelease = 'False'
      __version_release = 'False'

      >>> git log --oneline
      4e0596e Records tag for v1.2.0
      fc02d9d Bump version to v1.2.0

#. **ACT** to inform remote of bump

   * Temporarily disengage "Include administrators" on protected release branch.

    ::

      >>> git push origin release/1.2

      >>> git push origin v1.2.0

   * Now https://github.com/PCMSolver/pcmsolver/releases says ``v1.2.0`` & ``fc023d9d``

#. **EDIT** release description in the `GitHub web UI <https://github.com/PCMSolver/pcmsolver/releases>`_.

`Zenodo <https://zenodo.org/>`_ will automatically generate a new, versioned
DOI for the new release. It is no longer necessary to update the badge
in the ``README.md`` since it will always resolve to the latest released by
Zenodo.


How to create and remove an annotated Git tag on a remote
---------------------------------------------------------

PCMSolver versioning only works with *annotated* tags, not *lightweight*
tags as are created with the `GitHub interface
<https://github.com/PCMSolver/pcmsolver/releases/new>`_

* Create *annotated* tag::

    >>> git tag -a v1.1.12 <git hash if not current> -m "Version 1.1.12 released"
    >>> git push upstream --tags

* Delete tag::

    >>> git tag -d v1.1.12
    >>> git push origin :refs/tags/v1.1.12

* Pull tags::

    >>> git fetch <remote> 'refs/tags/*:refs/tags/*'
