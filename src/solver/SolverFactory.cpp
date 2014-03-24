#include "SolverFactory.hpp"

#include <stdexcept>
#include <string>

#include "Config.hpp"


class GreensFunction;

bool SolverFactory::registerSolver(std::string solverID,
                                   createSolverCallback createFunction)
{
    return callbacks.insert(CallbackMap::value_type(solverID, createFunction)).second;
}

bool SolverFactory::unRegisterSolver(std::string solverID)
{
    return callbacks.erase(solverID) == 1;
}

PCMSolver * SolverFactory::createSolver(std::string solverID,
                                        const solverData & _data)
{
    CallbackMap::const_iterator i = callbacks.find(solverID);
    if (i == callbacks.end()) {
        // The solverID was not found
        throw std::runtime_error("The unknown solver ID " + solverID +
                                 " occurred in SolverFactory.");
    }
    // Invoke the creation function
    return (i->second)(_data);
}
