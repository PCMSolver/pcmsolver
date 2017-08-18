/**
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

/*! \file Cxx11Workarounds.hpp
 *  \brief Provide hacks and workarounds for C++11
 *  \author Roberto Di Remigio
 *  \date 2015
 *
 *  As all hacks/workarounds this ain't a pretty piece of source code!
 */

/* --- Workarounds for constructs implemented by Boost prior to the C++11 standard
 * approval */
#ifdef HAS_CXX11
/* Smart pointers workarounds */
#include <memory>
namespace pcm {
using std::shared_ptr;
using std::make_shared;
using std::unique_ptr;
} /* end namespace pcm */
/* <functional> workarounds */
#include <functional>
namespace pcm {
using std::function;
using std::bind;
using std::placeholders::_1;
using std::placeholders::_2;
using std::placeholders::_3;
using std::placeholders::_4;
} /* end namespace pcm */
/* <tuple> workarounds */
#include <tuple>
namespace pcm {
using std::tuple;
using std::make_tuple;
using std::tie;
using std::ignore;
using std::get;
} /* end namespace pcm */
/* <array> workarounds */
#include <array>
namespace pcm {
using std::array;
} /* end namespace pcm */
/* std::to_string workarounds */
#include <string>
namespace pcm {
using std::to_string;
} /* end namespace pcm */
#else /* HAS_CXX11*/
/* Smart pointers workarounds */
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>
namespace pcm {
using boost::shared_ptr;
using boost::make_shared;
} /* end namespace pcm */
/* <functional> workarounds */
#include <boost/bind.hpp>
#include <boost/function.hpp>
namespace pcm {
using boost::function;
using boost::bind;
using ::_1;
using ::_2;
using ::_3;
using ::_4;
} /* end namespace pcm */
/* <tuple> workarounds */
#include <boost/tuple/tuple.hpp>
namespace pcm {
using boost::tuple;
using boost::make_tuple;
using boost::tie;
using boost::tuples::ignore;
using boost::get;
} /* end namespace pcm */
/* <array> workarounds */
#include <boost/array.hpp>
namespace pcm {
using boost::array;
} /* end namespace pcm */
/* std::to_string workarounds */
#include <boost/lexical_cast.hpp>
#include <string>
namespace pcm {
template <typename Source> std::string to_string(const Source & arg) {
  return boost::lexical_cast<std::string>(arg);
}
} /* end namespace pcm */
#endif /* HAS_CXX11 */

/* --- Workarounds for new keywords */
/* Workaroud for final */
#ifdef HAS_CXX11_CLASS_OVERRIDE
#define __final final
#else /* HAS_CXX11_CLASS_OVERRIDE */
#define __final
#endif /* HAS_CXX11_CLASS_OVERRIDE */

/* Workaroud for override */
#ifdef HAS_CXX11_CLASS_OVERRIDE
#define __override override
#else /* HAS_CXX11_CLASS_OVERRIDE */
#define __override
#endif /* HAS_CXX11_CLASS_OVERRIDE */

/* Workaroud for noexcept */
#ifdef HAS_CXX11_NOEXCEPT
#define __noexcept noexcept
#else /* HAS_CXX11_NOEXCEPT */
#define __noexcept throw()
#endif /* HAS_CXX11_NOEXCEPT */

/* Workaround for nullptr */
#ifdef HAS_CXX11_NULLPTR
#define __nullptr nullptr
#else /* HAS_CXX11_NULLPTR */
#define __nullptr NULL
#endif /* HAS_CXX11_NULLPTR */

/* Workaround for [[noreturn]] */
#ifdef HAS_CXX11_NORETURN
#define __noreturn [[noreturn]]
#elif defined(__GNUC__)
#define __noreturn __attribute__((noreturn))
#elif defined(_MSC_VER)
#define __declspec(noreturn)
#else /* HAS_CXX11_NORETURN */
#define __noreturn
#endif /* HAS_CXX11_NORETURN */
