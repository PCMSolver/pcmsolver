Solvers
=======

We will here describe the inheritance hierarchy for generating solvers, in order to use and extend it properly.
The runtime creation of solver objects relies on the Factory Method pattern [Gamma1994]_, [Alexandrescu2001]_, 
implemented through the SolverFactory class.

PCMSolver
---------
.. doxygenclass:: PCMSolver 
   :project: PCMSolver
   :members:
   :sections: public*, protected*, private*

IEFSolver
---------
.. doxygenclass:: IEFSolver 
   :project: PCMSolver
   :members:
   :sections: public*, protected*, private*

CPCMSolver
----------
.. doxygenclass:: CPCMSolver 
   :project: PCMSolver
   :members:
   :sections: public*, protected*, private*

WEMSolver
---------
.. doxygenclass:: WEMSolver 
   :project: PCMSolver
   :members:
   :sections: public*, protected*, private*

PWCSolver
---------
.. doxygenclass:: PWCSolver 
   :project: PCMSolver
   :members:
   :sections: public*, protected*, private*

PWLSolver
---------
.. doxygenclass:: PWLSolver 
   :project: PCMSolver
   :members:
   :sections: public*, protected*, private*

SolverFactory
-------------
.. doxygenclass:: SolverFactory
   :project: PCMSolver
   :members:
   :sections: public*, protected*, private*
