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

struct InputGePolTest {
    std::string filename;
    Input parsedInput;
    InputGePolTest() { SetUp(); }
    void SetUp() {
	filename = "gepol.inp";
	parsedInput = Input(filename);
    }
};

BOOST_FIXTURE_TEST_SUITE(InputGePol, InputGePolTest)

/*! \class Input 
 *  \test \b InputGePolTest_restart tests input reading by embedding Python pcmsolver.py script
 */
BOOST_FIXTURE_TEST_CASE(embdeddedPython, InputGePolTest)
{
	std::cout << "GePol" << std::endl;	
}

BOOST_AUTO_TEST_SUITE_END()

struct InputRestartTest {
    std::string filename;
    Input parsedInput;
    InputRestartTest() { SetUp(); }
    // List the contents of the input file here
    std::string cavFilename;
    void SetUp() {
	filename = "restart.inp";
	parsedInput = Input(filename);
        std::string cavFilename = "cavity.npz";
    }
};

BOOST_FIXTURE_TEST_SUITE(InputRestart, InputRestartTest)

/*! \class Input 
 *  \test \b InputRestartTest_restart tests input reading by embedding Python pcmsolver.py script
 */
BOOST_FIXTURE_TEST_CASE(embdeddedPython, InputRestartTest)
{
	std::cout << "Restart" << std::endl;	
}

BOOST_AUTO_TEST_SUITE_END()
