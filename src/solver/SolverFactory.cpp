#include <stdexcept>

#include "SolverFactory.h"

bool SolverFactory::registerSolver(std::string solverID, createSolverCallback createFunction)
{
	return callbacks.insert(CallbackMap::value_type(solverID, createFunction)).second;
}

bool SolverFactory::unRegisterSolver(std::string solverID)
{
	return callbacks.erase(solverID) == 1;
}

Solver * SolverFactory::createSolver(std::string solverID)
{
	CallbackMap::const_iterator i = callbacks.find(solverID);
	if (i == callbacks.end()) 
	{
		// The solverID was not found
                throw std::runtime_error("Unknown solver ID.");
	}
	// Invoke the creation function
	return (i->second)();
}
