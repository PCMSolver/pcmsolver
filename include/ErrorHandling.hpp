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

#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>

/*! \file ErrorHandling.hpp
 *  \brief Provide macros for error handling
 *  \author Roberto Di Remigio
 *  \date 2015
 *
 *  Exceptions:
 *
 *  \verbatim
 *  PCMSOLVER_ERROR(<Error Message>)
 *  \endverbatim
 *
 *  Use this to catch errors that might occur because of faulty
 *  data, i.e. other people's errors. Mainly in the API to the
 *  external world.
 *
 *  Assertions:
 *
 *  \verbatim
 *  PCMSOLVER_ASSERT(<Condition>)
 *  \endverbatim
 *
 *  This is just a wrapper around the standard assert macro.
 *  Use this to catch your own errors, i.e. broken preconditions
 *  to the internal functions/classes. In release mode assertions
 *  are not checked!
 *
 *  Static assertions:
 *
 *  \verbatim
 *  PCMSOLVER_STATIC_ASSERT(<Condition>, <Error Message>)
 *  \endverbatim
 *
 *  Uses static_assert. If not available, falls back to BOOST_STATIC_ASSERT_MSG
 *  Same usage pattern as for normal assertions. Static assertions are
 *  checked at compile-time.
 *  See also here:
 *http://www.boost.org/doc/libs/1_59_0/doc/html/boost_staticassert.html
 */

/// Macro to be used to signal error conditions
#define PCMSOLVER_ERROR(message)                                                    \
  {                                                                                 \
    std::ostringstream _err;                                                        \
    _err << "PCMSolver fatal error.\n"                                              \
         << " In function " << __func__ << " at line " << __LINE__ << " of file "   \
         << __FILE__ << "\n"                                                        \
         << message << std::endl;                                                   \
    std::fprintf(stderr, "%s\n", _err.str().c_str());                               \
    std::exit(EXIT_FAILURE);                                                        \
  }

/// Macro to be used for assertions
#define PCMSOLVER_ASSERT(arg) assert(arg)

/// Macro to be used for static assertions
#ifdef HAS_CXX11_STATIC_ASSERT
#define PCMSOLVER_STATIC_ASSERT(arg, msg) static_assert(arg, msg)
#else /* HAS_CXX11_STATIC_ASSERT */
#include <boost/static_assert.hpp>
#define PCMSOLVER_STATIC_ASSERT(arg, msg) BOOST_STATIC_ASSERT_MSG(arg, msg)
#endif /* HAS_CXX11_STATIC_ASSERT */
