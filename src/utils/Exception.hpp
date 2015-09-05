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

#ifndef EXCEPTION_HPP
#define EXCEPTION_HPP

#include <exception>
#include <string>

#include "Config.hpp"


/*! \file Exception.hpp
 *  \class Exception
 *  \brief Provide a simple exception class able to unwind the call stack
 *  \author Roberto Di Remigio
 *  \date 2015
 */

class Exception : public std::exception
{
    public:
        /*! Creates the exception object
         *  \param[in] m     error message
         *  \param[in] file  name of the file where the exception was thrown
         *  \param[in] line  number of the line where the exception was thrown
         *  \param[in] stack_size number of functions to be printed in the backtrace
         */
        Exception(const std::string & m, const char * file,
                const int line, const size_t stack_size = 5);
        virtual ~Exception() __noexcept {}
    private:
        /// The error message that will be printed out
        std::string message_;
        /// Function called when exception is thrown
        virtual const char *  what() const __noexcept { return message_.c_str(); }
};

/*! Returns call backtrace
 *  \param[out] size actual stack size
 *  \param[in] stack_size how far back in the call stack we want to go
 */
std::string unwind_call_stack(int size, const size_t stack_size = 5);

#endif // EXCEPTION_HPP
