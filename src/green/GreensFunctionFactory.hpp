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

#ifndef GREENSFUNCTIONFACTORY_HPP
#define GREENSFUNCTIONFACTORY_HPP

#include <string>
#include <map>

#include "Config.hpp"

#include "EigenPimpl.hpp"

struct greenData;
class IGreensFunction;

/*!
 *	\file GreensFunctionFactory.hpp
 *	\class GreensFunctionFactory
 *	\brief Implementation of the Factory Method for Green's functions.
 *	\author Roberto Di Remigio
 *	\date 2013
 *
 * 	Factory method implementation shamelessly copied from here \cite Alexandrescu2001
 * 	It is implemented as a Singleton.
 */

class GreensFunctionFactory
{
public:
    /*!
     * Callback function for Green's function creation.
     */
    typedef IGreensFunction * (*createGreensFunctionCallback)(const greenData & _data);
private:
    /*!
     * A map from the Green's function type identifier (a string) to its callback function.
     */
    typedef std::map<std::string, createGreensFunctionCallback> CallbackMap;
public:
    /*!
     * \brief Returns true if registration of the greenID was successful
     * \param greenID the Green's function identification string
     * \param createFunction the creation function related to the Green's function type given
     */
    bool registerGreensFunction(const std::string & greenID,
                                createGreensFunctionCallback createFunction);
    /*!
     * \brief Returns true if greenID was already registered
     * \param greenID the Green's function identification string
     */
    bool unRegisterGreensFunction(const std::string & greenID);
    /*!
     * Calls the appropriate creation function, based on the passed greenID
     */
    IGreensFunction * createGreensFunction(const std::string & greenID,
                                           const greenData & _data);
    /*!
     * Unique point of access to the unique instance of the GreensFunctionFactory
     */
    static GreensFunctionFactory& TheGreensFunctionFactory() {
        static GreensFunctionFactory obj;
        return obj;
    }
private:
    GreensFunctionFactory() {}
    /// Copy constructor is made private
    GreensFunctionFactory(const GreensFunctionFactory &other);
    GreensFunctionFactory& operator=(const GreensFunctionFactory &other);
    ~GreensFunctionFactory() {}
    CallbackMap callbacks;
};

#endif // GREENSFUNCTIONFACTORY_HPP
