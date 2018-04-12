Documentation
=============

This documentation is generated using `Sphinx <http://sphinx-doc.org/>`_ and
`Doxygen <http://www.stack.nl/~dimitri/doxygen/>`_ The two softwares are
bridged by means of the `Breathe extension <https://breathe.readthedocs.org/>`_
The online version of this documentation is built and served by `Read The Docs
<https://readthedocs.org/>`_.  The webpage http://pcmsolver.readthedocs.org/ is
updated on each push to the public GitHub repository.


How and what to document
------------------------

Doxygen enables documenting the code in the source code files thus removing a
"barrier" for developers.  To avoid that the code degenerates into a Big Ball
of Mud, it is mandatory to document directly within the source code classes and
functions. To document general programming principles, design choices,
maintenance etc. you can create a .rst file in the ``doc`` directory. Remember
to refer the new file inside the ``index.rst`` file (it won't be parsed
otherwise).  Sphing uses `reStructuredText
<http://docutils.sourceforge.net/rst.html>`_ and `Markdown
<https://daringfireball.net/projects/markdown/>`_. Support for Markdown is not
as extensive as for reStructuredText, see `these comments
<https://blog.readthedocs.com/adding-markdown-support/>`_. Follow the guidelines
in :cite:`Wilson2014` regarding what to document.

Write the documentation in the header file. To document a class, put
``/*! \class <myclass>`` inside the namespace but before the class.
Add the following to a ``.rst`` file:

.. code-block:: rst

  .. doxygenclass:: <namespace>::<myclass>
   :project: PCMSolver
   :members:
   :protected-members:
   :private-members:

Do similar when documenting ``struct``-s and complete files.

.. note::

   Use ``/*! */`` to open and close a Doxygen comment.

Documenting methods in derived classes
--------------------------------------

Virtual methods should only be documented in the base classes.
This avoids unnecessary verbosity and conforms to the principle: "Document
_what_, not _how_" :cite:`Wilson2014`
If you feel the _how_ needs to be explicitly documented, add some notes in the
appropriate ``.rst`` file.

How does this work?
-------------------

To have an offline version of the documentation just issue
in the ``doc`` folder:

.. code-block:: bash

   sphinx-build . _build

The HTML will be stored in ``_build/``. Open the ``_build/index.html`` file with
your browser to see and browse the documentation.

.. warning::

   It is only possible to build documentation locally from within the ``doc``
   folder.
   This choice was made to simplify the set up of the ReadTheDocs and local
   documentation building procedures and to minimize the chances of breaking
   either.

.. note::

   Building the documentation requires Python, Doxygen, Sphinx, Perl and the
   Python modules breathe, matplotlib, sphinx-rtd-theme, sphinxcontrib-bibtex
   and recommonmark.
   The required python modules can be installed by running ``pip install -r
   requirements.txt``.
   There is also a ``Pipfile`` in case people prefer to use ``pipenv``.
