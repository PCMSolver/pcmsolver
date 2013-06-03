#ifndef SOLVERFACTORY_H
#define SOLVERFACTORY_H

#include <iostream>
#include <string>
#include <map>

#include "Config.h"

#include <Eigen/Dense>


#include "PCMSolver.h"
#include "IEFSolver.h"
#include "CPCMSolver.h"
#include "WEMSolver.h"
#include "PWCSolver.h"
#include "PWLSolver.h"

/*
 * Factory method implementation shamelessly copied from
 * "Modern C++ Design" of A. Alexandrescu
 */


class SolverFactory {
	public:
		typedef PCMSolver * (*CreateSolverCallback)();
	private:
		typedef std::map<int, CreateSolverCallback> CallbackMap;
	public:
		// Returns 'true' if registration was successful
		bool RegisterSolver(int solverID, CreateSolverCallback createFunction);
		// Returns 'true' if registration it the solverType was registered before
		bool UnregisterSolver(int solverID);
		PCMSolver * CreateSolver(int solverID);
	private:
		CallbackMap callbacks;	
}

#endif
