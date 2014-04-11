Cavities
========

We will here describe the inheritance hierarchy for generating cavities, in order to use and extend it properly.
The runtime creation of cavity objects relies on the Factory Method pattern [Gamma1994]_, [Alexandrescu2001]_, 
implemented through the CavityFactory class.

Cavity
------
.. doxygenclass:: Cavity
   :project: PCMSolver
   :members:
   :sections: public*, protected*, private*

GePolCavity
-----------

.. doxygenclass:: GePolCavity 
   :project: PCMSolver
   :members:
   :sections: public*, protected*, private*

TsLessCavity
------------

.. doxygenclass:: TsLessCavity 
   :project: PCMSolver
   :members:
   :sections: public*, protected*, private*

RestartCavity
-------------

.. doxygenclass:: RestartCavity 
   :project: PCMSolver
   :members:
   :sections: public*, protected*, private*

WaveletCavity
-------------

.. doxygenclass:: WaveletCavity 
   :project: PCMSolver
   :members:
   :sections: public*, protected*, private*

CavityFactory
-------------

.. doxygenclass:: CavityFactory
   :project: PCMSolver
   :members:
   :sections: public*, protected*, private* 
