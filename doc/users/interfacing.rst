Interfacing a QM program and PCMSolver
======================================

.. warning::

   **Work in Progress!!**

How PCMSolver handles potentials and charges: surface functions
---------------------------------------------------------------

:cpp:class:`SurfaceFunction`

What you should care about: API functions
-----------------------------------------

These are the contents of the ``pcmsolver.h`` file defining
the public API of the PCMSolver library. The Fortran bindings
for the API are in the ``pcmsolver.F90`` file.
The indexing of symmetry operations and their mapping to a bitstring
is explained in the following Table. This is important when passing
symmetry information to the :cpp:func:`pcmsolver_new` function.

.. _symmetry-ops:
.. table:: Symmetry operations indexing within the module

   ===== === ========= ======
   Index zyx Generator Parity
   ===== === ========= ======
     0   000     E       1.0
     1   001    Oyz     -1.0
     2   010    Oxz     -1.0
     3   011    C2z      1.0
     4   100    Oxy     -1.0
     5   101    C2y      1.0
     6   110    C2x      1.0
     7   111     i      -1.0
   ===== === ========= ======


.. doxygenfile:: mock_pcmsolver.h
   :project: PCMSolver

Host input forwarding
---------------------

.. doxygenstruct:: PCMInput
   :project: PCMSolver

Internal details of the API
---------------------------

.. doxygenclass:: Meddle
   :project: PCMSolver
