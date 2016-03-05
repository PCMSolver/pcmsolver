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

#ifndef ERRORHANDLING_HPP
#define ERRORHANDLING_HPP

#include <cassert>
#include <cstdio>
#include <cstdlib>

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
 *  See also here: http://www.boost.org/doc/libs/1_59_0/doc/html/boost_staticassert.html
 */

/*! \brief Kills execution and prints out error message to stderr
 *  \param message Error message
 *  \param function Name of the function killing execution
 *  \param code Error code. Defaults to EXIT_FAILURE
 */
inline void pcmsolver_die(const std::string & message, const std::string & function, int code = EXIT_FAILURE)
{
  pcmsolver_die(message.c_str(), function.c_str(), code);
}

/*! \brief Kills execution and prints out error message to stderr
 *  \param message Error message
 *  \param function Name of the function killing execution
 *  \param code Error code. Defaults to EXIT_FAILURE
 */
inline void pcmsolver_die(const char * message, const char * function, int code = EXIT_FAILURE)
{
  std::fprintf(stderr, "In function: %s\n", function);
  std::fprintf(stderr, "PCMSolver fatal error %i: %s\n", code, message);
  std::exit(EXIT_FAILURE);
}

/// Macro to be used to signal error conditions
#define PCMSOLVER_ERROR(arg, func) pcmsolver_die(arg, func)

/// Macro to be used for assertions
#define PCMSOLVER_ASSERT(arg) assert(arg)

/// Macro to be used for static assertions
#ifdef HAS_CXX11_STATIC_ASSERT
#define PCMSOLVER_STATIC_ASSERT(arg, msg) static_assert(arg, msg)
#else /* HAS_CXX11_STATIC_ASSERT */
#include <boost/static_assert.hpp>
#define PCMSOLVER_STATIC_ASSERT(arg, msg) BOOST_STATIC_ASSERT_MSG(arg, msg)
#endif /* HAS_CXX11_STATIC_ASSERT */

#endif /* ERRORHANDLING_HPP */
