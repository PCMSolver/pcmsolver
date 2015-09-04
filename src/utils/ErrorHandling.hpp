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

#ifndef ERRORHANDLING_HPP
#define ERRORHANDLING_HPP

#include <cassert>

#include "Config.hpp"

#include <boost/static_assert.hpp>

#include "Exception.hpp"

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
 *  The output contains the most recent 5 function calls.
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
 *  Uses BOOST_STATIC_ASSERT_MSG, which falls back to the standard
 *  static_assert, if available.
 *  Same usage pattern as for normal assertions. Static assertions are
 *  checked at compile-time.
 *  See also here: http://www.boost.org/doc/libs/1_59_0/doc/html/boost_staticassert.html
 */

/// Macro to be used to throw exceptions
#define PCMSOLVER_ERROR(arg) throw Exception(arg, __FILE__, __LINE__)

/// Macro to be used for assertions
#define PCMSOLVER_ASSERT(arg) assert(arg)

/// Macro to be used for static assertions
#define PCMSOLVER_STATIC_ASSERT(arg, msg) BOOST_STATIC_ASSERT_MSG(arg, msg)

#endif // ERRORHANDLING_HPP
