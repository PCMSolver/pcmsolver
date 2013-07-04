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

/*!
 *	\file SolverFactory.h
 *	\class SolverFactory
 *
 *	\brief Implementation of the Factory Method for solvers. 
 *	\author Roberto Di Remigio
 *	\date 2013 
 *
 * 	Factory method implementation shamelessly copied from "Modern C++ Design" of A. Alexandrescu.
 * 	It is implemented as a Singleton.
 */

class SolverFactory {
	public:
		/*!
		 * Callback function for solver creation.
		 */
		typedef PCMSolver * (*createSolverCallback)();
	private:
		/*!
		 * A map from the solver type identifier (a string) to its callback function.
		 */
		typedef std::map<std::string, createSolverCallback> CallbackMap;
	public:
		/*!
		 * \brief Returns true if registration of the solverID was successful
		 * \param solverID the solver identification string
		 * \param createFunction the creation function related to the solver type given
		 */
		bool registerSolver(std::string solverID, createSolverCallback createFunction);
		/*!
		 * \brief Returns true if solverID was already registered
		 * \param solverID the solver identification string
		 */
		bool unRegisterSolver(std::string solverID);
		/*! 
		 * Calls the appropriate creation function, based on the passed cavityID
		 */
		PCMSolver * createSolver(std::string solverID, GreensFunction * gfInside_, GreensFunction * gfOutside_);
		/*!
		 * Unique point of access to the unique instance of the SolverFactory
		 */
		static SolverFactory& TheSolverFactory() 
		{
			static SolverFactory obj;
			return obj;
		}
	private:
		SolverFactory(){}
		/// Copy constructor is made private
		SolverFactory(const SolverFactory &other);
		SolverFactory& operator=(const SolverFactory &other);
        	~SolverFactory(){}
		CallbackMap callbacks;	
};

#endif
