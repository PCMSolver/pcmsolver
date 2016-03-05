/* pcmsolver_copyright_start */
/*
 *     PCMSolver, an API for the Polarizable Continuum Model
 *     Copyright (C) 2013-2015 Roberto Di Remigio, Luca Frediani and contributors
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
 *     PCMSolver API, see: <http://pcmsolver.readthedocs.org/>
 */
/* pcmsolver_copyright_end */

#ifndef FACTORY_HPP
#define FACTORY_HPP

#include <string>
#include <map>

#include "Config.hpp"

/*! \file Factory.hpp
 *	\class Factory
 *	\brief Implementation of the Factory Method
 *	\author Roberto Di Remigio
 *	\date 2015
 *  \tparam Object type of the object the factory will create
 *  \tparam ObjectInput type of the input wrapper struct
 *
 * 	Factory method implementation shamelessly copied from here \cite Alexandrescu2001
 * 	It is implemented as a Singleton.
 */

template <typename Object,
          typename ObjectInput>
class Factory __final
{
public:
    /*! \brief Callback function for object creation
     *  Returns a raw pointer of type Object; accepts an ObjectInput type
     */
    typedef pcm::function<Object * (const ObjectInput &)> creationalFunctor;
private:
    /*! std::map from the object type identifier (a string) to its callback function */
    typedef std::map<std::string, creationalFunctor> CallbackMap;
    /*! std::pair of an object type identifier and a callback function */
    typedef std::pair<std::string, creationalFunctor> CallbackPair;
public:
    /*! \brief Returns true on successful registration of the objID
     * \param[in] objID  the object's identification string
     * \param[in] functor the creation function related to the object type given
     */
    bool registerObject(const std::string & objID, const creationalFunctor & functor) {
        return callbacks_.insert(typename CallbackMap::value_type(objID, functor)).second;
    }
    /*! \brief Returns true if objID was already registered
     *  \param objID the object's identification string
     */
    bool unRegisterObject(const std::string & objID) {
        return callbacks_.erase(objID) == 1;
    }
    /*! \brief Calls the appropriate creation functor, based on the passed objID
     *  \param[in] objID the object's identification string
     *  \param[in] data  input data for the creation of the object
     */
    Object * create(const std::string & objID, const ObjectInput & data) {
        if (objID.empty()) PCMSOLVER_ERROR("No object identification string provided to the Factory.", BOOST_CURRENT_FUNCTION);
        typename CallbackMap::const_iterator i = callbacks_.find(objID);
        if (i == callbacks_.end()) PCMSOLVER_ERROR("The unknown object ID " + objID + " occurred in the Factory.", BOOST_CURRENT_FUNCTION);
        return (i->second)(data);
    }
    /*! Unique point of access to the unique instance of the Factory */
    static Factory & TheFactory() {
        static Factory obj;
        return obj;
    }
private:
    Factory() {}
    /// Copy constructor is made private
    Factory(const Factory & other);
    Factory& operator=(const Factory & other);
    ~Factory() {}
    CallbackMap callbacks_;
};

#endif // FACTORY_HPP
