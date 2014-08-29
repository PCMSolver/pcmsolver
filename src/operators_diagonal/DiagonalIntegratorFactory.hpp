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
 *     PCMSolver API, see: <http://pcmsolver.github.io/pcmsolver-doc>
 */
/* pcmsolver_copyright_end */

#ifndef DIAGONALINTEGRATORFACTORY_HPP
#define DIAGONALINTEGRATORFACTORY_HPP

#include <string>
#include <map>

#include "Config.hpp"

#include <Eigen/Dense>

class DiagonalIntegrator;

/*!
 *	\file DiagonalIntegratorFactory.hpp
 *	\class DiagonalIntegratorFactory
 *	\brief Implementation of the Factory Method for diagonal integrators 
 *	\author Roberto Di Remigio
 *	\date 2014
 *
 * 	Factory method implementation shamelessly copied from here \cite Alexandrescu2001
 * 	It is implemented as a Singleton.
 */

class DiagonalIntegratorFactory
{
public:
    /*!
     * Callback function for diagonal integrator creation.
     */
    typedef DiagonalIntegrator * (*createDiagonalIntegratorCallback)();
private:
    /*!
     * A map from the diagonal integrator type identifier (a string) to its callback function.
     */
    typedef std::map<std::string, createDiagonalIntegratorCallback> CallbackMap;
public:
    /*!
     * \brief Returns true if registration of the integratorID was successful
     * \param integratorID the diagonal integrator identification string
     * \param createFunction the creation function related to the diagonal integrator type given
     */
    bool registerDiagonalIntegrator(const std::string & integratorID,
                                createDiagonalIntegratorCallback createFunction);
    /*!
     * \brief Returns true if integratorID was already registered
     * \param integratorID the diagonal integrator identification string
     */
    bool unRegisterDiagonalIntegrator(const std::string & integratorID);
    /*!
     * Calls the appropriate creation function, based on the passed integratorID
     */
    DiagonalIntegrator * createDiagonalIntegrator(const std::string & integratorID);
    /*!
     * Unique point of access to the unique instance of the DiagonalIntegratorFactory
     */
    static DiagonalIntegratorFactory& TheDiagonalIntegratorFactory() {
        static DiagonalIntegratorFactory obj;
        return obj;
    }
private:
    DiagonalIntegratorFactory() {}
    /// Copy constructor is made private
    DiagonalIntegratorFactory(const DiagonalIntegratorFactory &other);
    DiagonalIntegratorFactory& operator=(const DiagonalIntegratorFactory &other);
    ~DiagonalIntegratorFactory() {}
    CallbackMap callbacks;
};

#endif // DIAGONALINTEGRATORFACTORY_HPP
