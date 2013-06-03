#include <stdexcept>

#include "SolverFactory.h"

;

bool SolverFactory::RegisterSolver(int solverID, CreateSolverCallback createFunction)
{
	return callbacks.insert(CallbackMap::value_type(solverID, createFunction)).second;
}

bool SolverFactory::UnregisterSolver(int solverID)
{
	return callbacks.erase(solverID) == 1;
}

PCMSolver * SolverFactory::CreateSolver(int solverID)
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
