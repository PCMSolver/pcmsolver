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

#ifndef CAVITYFACTORY_HPP
#define CAVITYFACTORY_HPP

#include <map>
#include <string>

#include "Config.hpp"


class Cavity;
struct cavityData;

/*!
 *	\file CavityFactory.hpp
 *	\class CavityFactory
 *	\brief Implementation of the Factory Method for cavities.
 *	\author Roberto Di Remigio
 *	\date 2013
 *
 * 	Factory method implementation shamelessly copied from here \cite Alexandrescu2001 
 * 	It is implemented as a Singleton.
 */

class CavityFactory
{
public:
    /*!
     * Callback function for cavity creation.
     */
    typedef Cavity * (*createCavityCallback)(const cavityData & _data);
private:
    /*!
     * A map from the cavity type identifier (a string) to its callback function.
     */
    typedef std::map<std::string, createCavityCallback> CallbackMap;
public:
    /*!
     * \brief Returns true if registration of the cavityID was successful
     * \param cavityID the cavity identification string
     * \param createFunction the creation function related to the cavity type given
     */
    bool registerCavity(std::string cavityID, createCavityCallback createFunction);
    /*!
     * \brief Returns true if cavityID was already registered
     * \param cavityID the cavity identification string
     */
    bool unRegisterCavity(std::string cavityID);
    /*!
     * Calls the appropriate creation function, based on the passed cavityID
     */
    Cavity * createCavity(std::string cavityID, const cavityData & _data);
    /*!
     * Unique point of access to the unique instance of the CavityFactory
     */
    static CavityFactory& TheCavityFactory() {
        static CavityFactory obj;
        return obj;
    }
private:
    CavityFactory() {}
    /// Copy constructor is made private
    CavityFactory(const CavityFactory &other);
    CavityFactory& operator=(const CavityFactory &other);
    ~CavityFactory() {}
    CallbackMap callbacks;
};

#endif // CAVITYFACTORY_HPP
