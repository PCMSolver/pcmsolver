Green's Functions
=================

We will here describe the inheritance hierarchy for generating Green's functions, in order to use and extend it properly.
The runtime creation of Green's functions objects relies on the Factory Method pattern [Gamma1994]_, [Alexandrescu2001]_, 
implemented through the GreensFunctionFactory class.

GreensFunctionInterface
-----------------------
.. doxygenclass:: GreensFunctionInterface
   :project: PCMSolver
   :members:
   :sections: public*, protected*, private*

GreensFunction
--------------
.. doxygenclass:: GreensFunction 
   :project: PCMSolver
   :members:
   :sections: public*, protected*, private*

GreensFunctionSum
-----------------
.. doxygenclass:: GreensFunctionSum
   :project: PCMSolver
   :members:
   :sections: public*, protected*, private*

Vacuum
------
.. doxygenclass:: Vacuum 
   :project: PCMSolver
   :members:
   :sections: public*, protected*, private*

UniformDielectric
-----------------
.. doxygenclass:: UniformDielectric 
   :project: PCMSolver
   :members:
   :sections: public*, protected*, private*

PlanarInterface
---------------
.. doxygenclass:: PlanarInterface 
   :project: PCMSolver
   :members:
   :sections: public*, protected*, private*

MetalSphere
-----------
.. doxygenclass:: MetalSphere 
   :project: PCMSolver
   :members:
   :sections: public*, protected*, private*

GreensFunctionFactory
---------------------
.. doxygenclass:: GreensFunctionFactory 
   :project: PCMSolver
   :members:
   :sections: public*, protected*, private*
