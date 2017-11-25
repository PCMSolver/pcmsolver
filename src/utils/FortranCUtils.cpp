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

#include "FortranCUtils.hpp"

#include <cctype>
#include <cstring>

#ifndef _fcdtocp
#define _fcdtocp(desc) (desc)
#endif

void pcmsolver_c2f_string(char * src, char * dest, int * len) {
  int sofar;
  for (sofar = 0; (sofar < *len) && (*src != '\0'); sofar++)
    *dest++ = *src++;
  while (sofar++ < *len)
    *dest++ = ' ';
}

void pcmsolver_f2c_string(char * src, char * dest, int * len) {
  char * str; /* Pointer to FORTRAN string */
  int i;      /* Local index variable */

  /* Search for the end of the string */
  str = _fcdtocp(src);
  for (i = *len - 1; i >= 0 && !std::isgraph((int)str[i]); i--)
    /*EMPTY*/;

  /* Copy text from FORTRAN to C string */
  std::memcpy(dest, str, (size_t)(i + 1));

  /* Terminate C string */
  dest[i + 1] = '\0';
}
