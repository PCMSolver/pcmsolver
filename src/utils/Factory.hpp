/*
 * PCMSolver, an API for the Polarizable Continuum Model
 * Copyright (C) 2017 Roberto Di Remigio, Luca Frediani and collaborators.
 *
 * This file is part of PCMSolver.
 *
 * PCMSolver is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * PCMSolver is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with PCMSolver.  If not, see <http://www.gnu.org/licenses/>.
 *
 * For information on the complete list of contributors to the
 * PCMSolver API, see: <http://pcmsolver.readthedocs.io/>
 */

#pragma once

#include <map>
#include <string>

#include "Config.hpp"

#ifdef HAS_CXX11_VARIADIC_TEMPLATES
#include <type_traits>
#else /* HAS_CXX11_VARIADIC_TEMPLATES */
#include <boost/utility/result_of.hpp>
#endif /* HAS_CXX11_VARIADIC_TEMPLATES */

namespace pcm {
namespace detail {
/*! \file Factory.hpp
 *  \class BaseFactory
 *  \brief A base class for the Factory Method
 *  \author Roberto Di Remigio
 *  \date 2017
 *  \tparam CreateObject type of the object creation callback function
 *
 *  This class provides basic functionality to implement the Factory Method
 *  as described by Alexandrescu in \cite Alexandrescu2001
 */
template <typename CreateObject> class BaseFactory {
private:
  /*! std::map from the object type identifier (a string) to its callback function */
  typedef std::map<std::string, CreateObject> CallbackMap;
  /*! std::pair of an object type identifier and a callback function */
  typedef typename CallbackMap::value_type CallbackPair;
  /*! const iterator */
  typedef typename CallbackMap::const_iterator CallbackConstIter;

protected:
  CallbackMap callbacks_;
  /*! \brief Retrieve constant iterator from map given object identifier
   * \param[in] objID  the object's identification string
   */
  CallbackConstIter retrieve(const std::string & objID) const {
    if (objID.empty())
      PCMSOLVER_ERROR("No object identification string provided to the Factory.");
    CallbackConstIter i = callbacks_.find(objID);
    if (i == callbacks_.end())
      PCMSOLVER_ERROR("The unknown object ID " + objID +
                      " occurred in the Factory.");
    return i;
  }

private:
  /*! \brief Returns true on successful registration of the objID
   * \param[in] objID  the object's identification string
   * \param[in] functor the creation function related to the object type given
   */
  bool registerObject(const std::string & objID, const CreateObject & functor) {
    return callbacks_.insert(CallbackPair(objID, functor)).second;
  }
  /*! \brief Returns true if objID was already registered
   *  \param objID the object's identification string
   */
  bool unRegisterObject(const std::string & objID) {
    return callbacks_.erase(objID) == 1;
  }

public:
  /*! \brief Subscribes object with objID to factory
   * \param[in] objID  the object's identification string
   * \param[in] functor the creation function related to the object type given
   */
  void subscribe(const std::string & objID, const CreateObject & functor) {
    bool done = this->registerObject(objID, functor);
    if (!done)
      PCMSOLVER_ERROR("Subscription of object ID " + objID + " to factory failed!");
  }
  /*! \brief Unsubscribes object with objID from factory
   *  \param objID the object's identification string
   */
  void unsubscribe(const std::string & objID) {
    bool done = this->unRegisterObject(objID);
    if (!done)
      PCMSOLVER_ERROR("Unsubscription of object ID " + objID +
                      " from factory failed!");
  }
};
} // namespace detail
#ifdef HAS_CXX11_VARIADIC_TEMPLATES
/*! \file Factory.hpp
 *  \class Factory
 *  \brief C++11 implementation of the Factory Method
 *  \author Roberto Di Remigio
 *  \date 2017
 *  \tparam CreateObject type of the object creation callback function
 *
 *  Factory Method using a variadic template parameter pack to handle possible
 *  input value to the object creation callback function.
 */
template <typename CreateObject>
class Factory __final : public detail::BaseFactory<CreateObject> {
public:
  /*! \brief Calls the appropriate creation functor, based on the passed objID
   *  \param[in] objID the object's identification string
   *  \param[in] data  input data for the creation of the object
   *  \tparam ObjectInputArgs type(s) of the object creation callback function input
   * arguments
   *  \note Return type is deduced based on the type(s) of the input
   *  argument(s) template parameter pack and the type of the object creation
   *  callback function.
   */
  template <typename... ObjectInputArgs>
  typename std::result_of<CreateObject(ObjectInputArgs...)>::type create(
      const std::string & objID,
      ObjectInputArgs... data) const {
    return (this->retrieve(objID)->second)(data...);
  }
};
#else  /* HAS_CXX11_VARIADIC_TEMPLATES */
/*! \file Factory.hpp
 *  \class Factory
 *  \brief C++03 implementation of the one- or zero-argument Factory Method
 *  \author Roberto Di Remigio
 *  \date 2015-2017
 *  \tparam CreateObject type of the object creation callback function
 *
 *  \warning This will only work when the CreateObject function accepts zero or
 *  one input arguments.
 */
template <typename CreateObject>
class Factory : public detail::BaseFactory<CreateObject> {
public:
  /*! \brief Calls the appropriate creation functor, based on the passed objID
   *  \param[in] objID the object's identification string
   *  \param[in] data  input data for the creation of the object
   *  \note This is the one-parameter version. Return type is deduced based on
   *  the type of the input argument template parameter and the type of the
   *  object creation callback function.
   */
  template <typename ObjectInputArg>
  typename boost::result_of<CreateObject(ObjectInputArg)>::type create(
      const std::string & objID,
      const ObjectInputArg & data) const {
    return (this->retrieve(objID)->second)(data);
  }

  /*! \brief Calls the appropriate creation functor, based on the passed objID
   *  \param[in] objID the object's identification string
   *  \note This is the zero-parameter version. Return type is deduced based on
   *  the type of the type of the object creation callback function.
   */
  typename boost::result_of<CreateObject()>::type create(
      const std::string & objID) const {
    return (this->retrieve(objID)->second)();
  }
};
#endif /* HAS_CXX11_VARIADIC_TEMPLATES */
} // namespace pcm
