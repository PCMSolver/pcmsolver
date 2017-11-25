Interfacing a QM program and PCMSolver
======================================

For the impatients: tl;dr
-------------------------

In these examples, we want to show how *every function* in the API works.
If your program is written in Fortran, head over to :ref:`fortran-example`
If your program is written in C/C++, head over to :ref:`C-example`

How PCMSolver handles potentials and charges: surface functions
---------------------------------------------------------------

Electrostatic potential vectors and the corresponding apparent surface
charge vectors are handled internally as `surface functions`.
The actual values are stored into Eigen vectors and saved into a
map. The mapping is between the name of the surface function, given by
the programmer writing the interface to the library, and the vector holding
the values.

What you should care about: API functions
-----------------------------------------

These are the contents of the ``pcmsolver.h`` file defining
the public API of the PCMSolver library. The Fortran bindings
for the API are in the ``pcmsolver.f90`` file.
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

.. doxygenclass:: pcm::Meddle
   :project: PCMSolver
   :members:
   :protected-members:
   :private-members:

.. doxygenclass:: pcm::Input
   :project: PCMSolver
   :members:
   :protected-members:
   :private-members:
