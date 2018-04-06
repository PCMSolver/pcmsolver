Green's Functions
=================

We will here describe the inheritance hierarchy for generating Green's
functions, in order to use and extend it properly.  The runtime creation of
Green's functions objects relies on the Factory Method pattern
:cite:`Gamma1994,Alexandrescu2001`, implemented through the
generic Factory class.

The top-level header, _i.e._ to be included in client code, is ``Green.hpp``.
The common interface to all Green's function classes is specified by the ``IGreensFunction`` class,
this is non-templated.
All other classes are templated.
The Green's functions are registered to the factory based on a label encoding: type, derivative, and dielectric profile.
The only allowed labels must be listed in ``src/green/Green.hpp``. If they are not, they can not be selected at run time.

.. image:: ../gfx/bar_charts/green.svg
   :scale: 70 %
   :align: center

IGreensFunction
---------------

.. doxygenclass:: pcm::IGreensFunction
   :project: PCMSolver
   :members:
   :protected-members:
   :private-members:

GreensFunction
--------------
.. doxygenclass:: pcm::green::GreensFunction
   :project: PCMSolver
   :members:
   :protected-members:
   :private-members:

Vacuum
------
.. doxygenclass:: pcm::green::Vacuum
   :project: PCMSolver
   :members:
   :protected-members:
   :private-members:

UniformDielectric
-----------------
.. doxygenclass:: pcm::green::UniformDielectric
   :project: PCMSolver
   :members:
   :protected-members:
   :private-members:

IonicLiquid
-----------
.. doxygenclass:: pcm::green::IonicLiquid
   :project: PCMSolver
   :members:
   :protected-members:
   :private-members:

AnisotropicLiquid
-----------------
.. doxygenclass:: pcm::green::AnisotropicLiquid
   :project: PCMSolver
   :members:
   :protected-members:
   :private-members:

SphericalDiffuse
----------------
.. doxygenclass:: pcm::green::SphericalDiffuse
   :project: PCMSolver
   :members:
   :protected-members:
   :private-members:
