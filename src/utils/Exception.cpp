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

#include "Exception.hpp"

#include <cstdlib>
#include <sstream>
#include <string>
#include <vector>
#include <execinfo.h>
#include <cxxabi.h>

#include "Config.hpp"

inline std::string unwind_call_stack(int size, const size_t stack_size) {
    std::stringstream info;
    std::vector<void *> Stack(stack_size);
    char ** strings;
    size = backtrace(&Stack[0], stack_size);
    int status = -1;
    strings = backtrace_symbols(&Stack[0], size);
    for (int i = 0; i < size; ++i) {
      //This part from https://panthema.net/2008/0901-stacktrace-demangled/
      char *begin_name = NULL, *begin_offset = NULL, *end_offset =NULL;
      for (char *p = strings[i]; *p; ++p) {
         if (*p == '(') begin_name = p;
         else if (*p == '+') begin_offset = p;
         else if (*p == ')' && begin_offset) {
            end_offset = p;
            break;
         }
      }
      if (begin_name && begin_offset && end_offset && begin_name<begin_offset) {
             *begin_name++ = '\0';
             *begin_offset++ = '\0';
             *end_offset = '\0';
             char * demangled = abi::__cxa_demangle(begin_name, 0, 0, &status);
             if (status == 0) info << demangled << std::endl;
             free(demangled);
      }
    }
    return info.str();
}

Exception::Exception(const std::string & m, const char * file, const int line, const size_t stack_size) {
    std::stringstream error;
    error << std::endl << "Fatal error: " << m << std::endl;
    error << "Error occurred in file: " << file << " on line: " << line << std::endl;

    size_t size = 0;
    std::string call_stack_info = unwind_call_stack(size, stack_size);

    error << "The most recent " << (size < stack_size ? size : stack_size)
          << " function calls were:" << std::endl << std::endl;

    message_ = error.str() + call_stack_info;
}
