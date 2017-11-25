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

#include "FCMangle.hpp"

#ifdef __cplusplus
extern "C" {
#endif

/*! \brief Convert a C string to a Fortran string
 *  \param[in] src the string to convert
 *  \param[out] dest the converted string
 *  \param[in]      len length of Fortran string
 *  Change a C string (NULL terminated) into a Fortran string.
 *  Basically, all that is done is that the NULL is ripped out
 *  and the string is padded with spaces
 */
#define pcmsolver_c2f_string                                                        \
  FortranCInterface_GLOBAL_(pcmsolver_c2f_string, PCMSOLVER_C2F_STRING)
void pcmsolver_c2f_string(char * src, char * dest, int * len);

/*! \brief Convert a Fortran string to a C string
 *  \param[in] src Fortran string
 *  \param[out] dest C string
 *  \param[in] len length of Fortran string
 *  Chop off trailing blanks off of a Fortran string and
 *  move it into C string.
 */
#define pcmsolver_f2c_string                                                        \
  FortranCInterface_GLOBAL_(pcmsolver_f2c_string, PCMSOLVER_F2C_STRING)
void pcmsolver_f2c_string(char * src, char * dest, int * len);

#ifdef __cplusplus
}
#endif
