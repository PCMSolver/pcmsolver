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

#define BOOST_TEST_MODULE Input 

#include <iostream>
#include <string>

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include "Config.hpp"

#include <Eigen/Dense>

#include "Input.hpp"

struct InputTest {
    std::string filename;
    Input * parsedInput;
    InputTest() { SetUp(); }
    void SetUp() {
	filename = "parsingTest.inp";
	parsedInput = &Input::TheInput(filename);
    }
};

/*! \class Input 
 *  \test \b InputTest_embeddedPython tests input reading using Python embedding
 */
BOOST_FIXTURE_TEST_CASE(embdeddedPython, InputTest)
{
	std::cout << "BLA BLA" << std::endl;	
}
