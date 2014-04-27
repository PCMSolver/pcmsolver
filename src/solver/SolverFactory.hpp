/* pcmsolver_copyright_start */
/*
 *     PCMSolver, an API for the Polarizable Continuum Model
 *     Copyright (C) 2013 Roberto Di Remigio, Luca Frediani and contributors
 *     
 *     This file is part of PCMSolver.
 *
 *     PCMSolver is free software: you can redistribute it and/or modify       
 *     it under the terms of the GNU Lesser General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *                                                                          
 *     PCMSolver is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU Lesser General Public License for more details.
 *                                                                          
 *     You should have received a copy of the GNU Lesser General Public License
 *     along with PCMSolver.  If not, see <http://www.gnu.org/licenses/>.
 *
 *     For information on the complete list of contributors to the
 *     PCMSolver API, see: <https://repo.ctcc.no/projects/pcmsolver>
 */
/* pcmsolver_copyright_end */

#ifndef SOLVERFACTORY_HPP
#define SOLVERFACTORY_HPP

#include <string>
#include <map>

#include "Config.hpp"


struct solverData;
class PCMSolver;

/*!
 *	\file SolverFactory.hpp
 *	\class SolverFactory
 *	\brief Implementation of the Factory Method for solvers.
 *	\author Roberto Di Remigio
 *	\date 2013
 *
 * 	Factory method implementation shamelessly copied from here \cite Alexandrescu2001 
 * 	It is implemented as a Singleton.
 */

class SolverFactory
{
public:
    /*!
     * Callback function for solver creation.
     */
    typedef PCMSolver * (*createSolverCallback)(const solverData & _data);
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
    PCMSolver * createSolver(std::string solverID, const solverData & _data);
    /*!
     * Unique point of access to the unique instance of the SolverFactory
     */
    static SolverFactory& TheSolverFactory() {
        static SolverFactory obj;
        return obj;
    }
private:
    SolverFactory() {}
    /// Copy constructor is made private
    SolverFactory(const SolverFactory &other);
    SolverFactory& operator=(const SolverFactory &other);
    ~SolverFactory() {}
    CallbackMap callbacks;
};

#endif // SOLVERFACTORY_HPP
