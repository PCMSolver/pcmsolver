/*
 * PCMSolver, an API for the Polarizable Continuum Model
 * Copyright (C) 2019 Roberto Di Remigio, Luca Frediani and contributors.
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
using std::make_shared;
using std::shared_ptr;
using std::unique_ptr;
} /* end namespace pcm */
/* <functional> workarounds */
#include <functional>
namespace pcm {
using std::bind;
using std::function;
using std::placeholders::_1;
using std::placeholders::_2;
using std::placeholders::_3;
using std::placeholders::_4;
} /* end namespace pcm */
/* <tuple> workarounds */
#include <tuple>
namespace pcm {
using std::get;
using std::ignore;
using std::make_tuple;
using std::tie;
using std::tuple;
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
/* <cmath> workarounds */
#include <cmath>
namespace pcm {
using std::erf;
} /* end namespace pcm */
/* <type_traits> workarounds */
#include <type_traits>
namespace pcm {
using std::enable_if;
using std::is_same;
} /* end namespace pcm */
#else /* HAS_CXX11*/
/* Smart pointers workarounds */
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>
namespace pcm {
using boost::make_shared;
using boost::shared_ptr;
} /* end namespace pcm */
/* <functional> workarounds */
#include <boost/bind.hpp>
#include <boost/function.hpp>
namespace pcm {
using ::_1;
using ::_2;
using ::_3;
using ::_4;
using boost::bind;
using boost::function;
} /* end namespace pcm */
/* <tuple> workarounds */
#include <boost/tuple/tuple.hpp>
namespace pcm {
using boost::get;
using boost::make_tuple;
using boost::tie;
using boost::tuple;
using boost::tuples::ignore;
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
/* <cmath> workarounds */
#include <boost/math/special_functions/erf.hpp>
namespace pcm {
using boost::math::erf;
} /* end namespace pcm */
/* <type_traits> workarounds */
#include <boost/core/enable_if.hpp>
#include <boost/type_traits/is_same.hpp>
namespace pcm {
using boost::enable_if;
using boost::is_same;
} /* end namespace pcm */
#endif /* HAS_CXX11 */

/* --- Workarounds for new keywords */
/* Workaroud for final */
#ifdef HAS_CXX11
#define __final final
#else /* HAS_CXX11 */
#define __final
#endif /* HAS_CXX11 */

/* Workaroud for override */
#ifdef HAS_CXX11
#define __override override
#else /* HAS_CXX11 */
#define __override
#endif /* HAS_CXX11 */

/* Workaroud for noexcept */
#ifdef HAS_CXX11
#define __noexcept noexcept
#else /* HAS_CXX11 */
#define __noexcept throw()
#endif /* HAS_CXX11 */

/* Workaround for nullptr */
#ifdef HAS_CXX11
#define __nullptr nullptr
#else /* HAS_CXX11 */
#define __nullptr NULL
#endif /* HAS_CXX11 */

/* Workaround for [[noreturn]] */
#ifdef HAS_CXX11
#define __noreturn [[noreturn]]
#elif defined(__GNUC__)
#define __noreturn __attribute__((noreturn))
#elif defined(_MSC_VER)
#define __declspec(noreturn)
#else /* HAS_CXX11 */
#define __noreturn
#endif /* HAS_CXX11 */
