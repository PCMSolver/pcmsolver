/* pcmsolver_copyright_start */
/*
 *     PCMSolver, an API for the Polarizable Continuum Model
 *     Copyright (C) 2013-2016 Roberto Di Remigio, Luca Frediani and contributors
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
 *     PCMSolver API, see: <http://pcmsolver.readthedocs.io/>
 */
/* pcmsolver_copyright_end */

#define CATCH_CONFIG_RUNNER
#include "catch.hpp"

#include "PhysicalConstants.hpp"

extern "C"
void host_writer(const char * /* message */, int /* message_length */);

int main( int argc, char* const argv[] )
{
  // global setup...
  initBohrToAngstrom(bohrToAngstrom);

  int result = Catch::Session().run( argc, argv );

  // global clean-up...

  return result;
}

extern "C"
void host_writer(const char * /* message */, int /* message_length */) {}
