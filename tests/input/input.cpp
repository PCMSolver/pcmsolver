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
#include "PhysicalConstants.hpp"

struct InputGePolTest {
    std::string filename;
    Input parsedInput;
    InputGePolTest() { SetUp(); }
    // List the contents of the input file here
    int CODATAyear;
    std::string type;
    std::string cavFilename;
    int patchLevel;
    double coarsity;
    double area;
    double minDistance;
    int derOrder;
    bool scaling;
    std::string radiiSet;
    double minimalRadius;
    std::string mode;
    Solvent solvent;
    bool hasSolvent;
    std::string solverType;
    int equationType;
    double correction;
    bool hermitivitize;
    double probeRadius;
    std::string greenInsideType;
    std::string greenOutsideType;
    int derivativeInsideType;
    int derivativeOutsideType;
    double epsilonInside;
    double epsilonOutside;
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
    int CODATAyear;
    std::string type;
    std::string cavFilename;
    int patchLevel;
    double coarsity;
    double area;
    double minDistance;
    int derOrder;
    bool scaling;
    std::string radiiSet;
    double minimalRadius;
    std::string mode;
    std::vector<int> atoms;
    std::vector<double> radii;
    std::vector<Sphere> spheres;
    std::string solvent;
    std::string solverType;
    int equationType;
    double correction;
    bool hermitivitize;
    double probeRadius;
    std::string greenInsideType;
    std::string greenOutsideType;
    int derivativeInsideType;
    int derivativeOutsideType;
    double epsilonInside;
    double epsilonOutside;
    void SetUp() {
	filename = "restart.inp";
	parsedInput = Input(filename);
	CODATAyear = 2010;
        type = "RESTART";
        cavFilename = "cavity.npz";
	patchLevel  = 2;
	coarsity = 0.5;
	area = 0.3;
	minDistance = 0.1;
	derOrder = 4;
        scaling = true;
	radiiSet = "BONDI";
	minimalRadius = 100.0;
	mode = "IMPLICIT";
	solvent = "Water"; // Name in the Solvent object
	solverType = "IEFPCM";
	equationType = 1; 
	correction = 0.0;
	hermitivitize = true;
	probeRadius = 1.385 * angstromToBohr(CODATAyear); // The value for water
	greenInsideType = "VACUUM";
	greenOutsideType = "UNIFORMDIELECTRIC";
	derivativeInsideType = 1;
	derivativeOutsideType = 1;
    }
};

BOOST_FIXTURE_TEST_SUITE(InputRestart, InputRestartTest)

/*! \class Input 
 *  \test \b InputRestartTest_restart tests input reading by embedding Python pcmsolver.py script
 */
BOOST_FIXTURE_TEST_CASE(embdeddedPython, InputRestartTest)
{
    double threshold = 1.0e-12;
    BOOST_REQUIRE_EQUAL(CODATAyear,            parsedInput.CODATAyear());           
    BOOST_REQUIRE_EQUAL(type,                  parsedInput.cavityType());
    BOOST_REQUIRE_EQUAL(cavFilename,           parsedInput.cavityFilename());
    BOOST_REQUIRE_EQUAL(patchLevel,            parsedInput.patchLevel());
    BOOST_REQUIRE_CLOSE(coarsity,              parsedInput.coarsity(), threshold);
    BOOST_REQUIRE_CLOSE(area,                  parsedInput.area(), threshold);
    BOOST_REQUIRE_CLOSE(minDistance,           parsedInput.minDistance(), threshold);
    BOOST_REQUIRE_EQUAL(derOrder,              parsedInput.derOrder());
    BOOST_REQUIRE_EQUAL(scaling,               parsedInput.scaling());
    BOOST_REQUIRE_EQUAL(radiiSet,              parsedInput.radiiSet());
    BOOST_REQUIRE_CLOSE(minimalRadius,         parsedInput.minimalRadius(), threshold);
    BOOST_REQUIRE_EQUAL(mode,                  parsedInput.mode());
    BOOST_REQUIRE_EQUAL(solvent,               parsedInput.solvent().name());
    BOOST_REQUIRE_EQUAL(solverType,            parsedInput.solverType());
    BOOST_REQUIRE_EQUAL(equationType,          parsedInput.equationType());
    BOOST_REQUIRE_CLOSE(correction,            parsedInput.correction(), threshold);
    BOOST_REQUIRE_EQUAL(hermitivitize,         parsedInput.hermitivitize());
    BOOST_REQUIRE_CLOSE(probeRadius,           parsedInput.probeRadius(), threshold);
    BOOST_REQUIRE_EQUAL(greenInsideType,       parsedInput.greenInsideType());
    BOOST_REQUIRE_EQUAL(greenOutsideType,      parsedInput.greenOutsideType());
    BOOST_REQUIRE_EQUAL(derivativeInsideType,  parsedInput.derivativeInsideType());
    BOOST_REQUIRE_EQUAL(derivativeOutsideType, parsedInput.derivativeOutsideType());
}

BOOST_AUTO_TEST_SUITE_END()
