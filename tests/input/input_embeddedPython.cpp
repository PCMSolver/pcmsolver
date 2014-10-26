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
#include <vector>

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include "Config.hpp"

#include <Eigen/Dense>

#include "Input.hpp"
#include "PhysicalConstants.hpp"
#include "Sphere.hpp"

struct InputTsLessTest {
    std::string filename;
    Input parsedInput;
    InputTsLessTest() { SetUp(); }
    // List the contents of the input file here
    std::string units;
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
	filename = "tsless.inp";
	parsedInput = Input(filename);
	units = "ANGSTROM";
	CODATAyear = 2002;
        type = "TSLESS";
	area = 0.6 * angstrom2ToBohr2(CODATAyear);
	minDistance = 5.0 * angstromToBohr(CODATAyear);
	derOrder = 25;
        scaling = false;
	radiiSet = "UFF";
	minimalRadius = 0.1 * angstromToBohr(CODATAyear);
	mode = "ATOMS";
	atoms.push_back(1); atoms.push_back(2); atoms.push_back(4);
	radii.push_back(1.4 * angstromToBohr(CODATAyear));
	radii.push_back(4.3 * angstromToBohr(CODATAyear));
	radii.push_back(1.2 * angstromToBohr(CODATAyear));
	solvent = "Cyclohexane"; // Name in the Solvent object
	solverType = "CPCM";
	correction = 0.5;
	hermitivitize = false;
	probeRadius = 2.815 * angstromToBohr(CODATAyear); // The value for water
	greenInsideType = "VACUUM";
	greenOutsideType = "UNIFORMDIELECTRIC";
	derivativeInsideType = 1;
	derivativeOutsideType = 1;
    }
};

BOOST_FIXTURE_TEST_SUITE(InputTsLess, InputTsLessTest)

/*! \class Input 
 *  \test \b InputTsLessTest_embeddedPythonTsLess tests input reading by embedding Python pcmsolver.py script
 */
BOOST_FIXTURE_TEST_CASE(embeddedPythonTsLess, InputTsLessTest)
{
    BOOST_REQUIRE_EQUAL(units,                 parsedInput.units());
    BOOST_REQUIRE_EQUAL(CODATAyear,            parsedInput.CODATAyear());           
    BOOST_REQUIRE_EQUAL(type,                  parsedInput.cavityType());
    BOOST_REQUIRE_CLOSE(area,                  parsedInput.area(), 1.0e-11);
    BOOST_REQUIRE_CLOSE(minDistance,           parsedInput.minDistance(), 1.0e-10);
    BOOST_REQUIRE_EQUAL(derOrder,              parsedInput.derOrder());
    BOOST_REQUIRE_EQUAL(scaling,               parsedInput.scaling());
    BOOST_REQUIRE_EQUAL(radiiSet,              parsedInput.radiiSet());
    BOOST_REQUIRE_CLOSE(minimalRadius,         parsedInput.minimalRadius(), 5.0e-10);
    BOOST_REQUIRE_EQUAL(mode,                  parsedInput.mode());
    for (size_t i = 0; i < atoms.size(); ++i) {
	    BOOST_REQUIRE_EQUAL(atoms[i],      parsedInput.atoms(i));
	    BOOST_REQUIRE_CLOSE(radii[i],      parsedInput.radii(i), 5.0e-10);
    }
    BOOST_REQUIRE_EQUAL(solvent,               parsedInput.solvent().name());
    BOOST_REQUIRE_EQUAL(solverType,            parsedInput.solverType());
    BOOST_REQUIRE_CLOSE(correction,            parsedInput.correction(), 1.0e-12);
    BOOST_REQUIRE_EQUAL(hermitivitize,         parsedInput.hermitivitize());
    BOOST_REQUIRE_CLOSE(probeRadius,           parsedInput.probeRadius(), 1.0e-12);
    BOOST_REQUIRE_EQUAL(greenInsideType,       parsedInput.greenInsideType());
    BOOST_REQUIRE_EQUAL(greenOutsideType,      parsedInput.greenOutsideType());
    BOOST_REQUIRE_EQUAL(derivativeInsideType,  parsedInput.derivativeInsideType());
    BOOST_REQUIRE_EQUAL(derivativeOutsideType, parsedInput.derivativeOutsideType());
}

BOOST_AUTO_TEST_SUITE_END()

struct InputRestartTest {
    std::string filename;
    Input parsedInput;
    InputRestartTest() { SetUp(); }
    // List the contents of the input file here
    std::string units;
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
	units = "AU";
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
 *  \test \b InputRestartTest_embeddedPythonRestart tests input reading by embedding Python pcmsolver.py script
 */
BOOST_FIXTURE_TEST_CASE(embeddedPythonRestart, InputRestartTest)
{
    double threshold = 1.0e-12;
    BOOST_REQUIRE_EQUAL(units,                 parsedInput.units());
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

struct InputWaveletTest {
    std::string filename;
    Input parsedInput;
    InputWaveletTest() { SetUp(); }
    // List the contents of the input file here
    std::string units;
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
	filename = "wavelet.inp";
	parsedInput = Input(filename);
	units = "ANGSTROM";
	CODATAyear = 1998;
        type = "WAVELET";
	patchLevel = 1;
	coarsity   = 0.3;
        scaling = true;
	radiiSet = "BONDI";
	mode = "EXPLICIT";
	Eigen::Vector3d c1, c2, c3;
	c1 <<  0.00,  0.00, 0.00;
	c2 <<  0.00, -0.96, 0.00;
	c3 << -0.905, 0.32, 0.00;
	c2 *= angstromToBohr(CODATAyear);
	c3 *= angstromToBohr(CODATAyear);
	Sphere sph1(c1, 1.80 * angstromToBohr(CODATAyear));
	Sphere sph2(c2, 1.44 * angstromToBohr(CODATAyear));
	Sphere sph3(c3, 1.44 * angstromToBohr(CODATAyear));
	solverType = "WAVELET";
	equationType = 0;
	probeRadius = 1.385 * angstromToBohr(CODATAyear); // The value for water
	greenInsideType = "VACUUM";
	greenOutsideType = "UNIFORMDIELECTRIC";
	derivativeInsideType = 0;
	derivativeOutsideType = 2;
	epsilonInside = 1.0;
	epsilonOutside = 78.39;
    }
};

BOOST_FIXTURE_TEST_SUITE(InputWavelet, InputWaveletTest)

/*! \class Input 
 *  \test \b InputWaveletTest_embeddedPythonWavelet tests input reading by embedding Python pcmsolver.py script
 */
BOOST_FIXTURE_TEST_CASE(embeddedPythonWavelet, InputWaveletTest)
{
    double threshold = 1.0e-12;
    BOOST_REQUIRE_EQUAL(units,                 parsedInput.units());
    BOOST_REQUIRE_EQUAL(CODATAyear,            parsedInput.CODATAyear());           
    BOOST_REQUIRE_EQUAL(type,                  parsedInput.cavityType());
    BOOST_REQUIRE_EQUAL(patchLevel,            parsedInput.patchLevel());
    BOOST_REQUIRE_CLOSE(coarsity,              parsedInput.coarsity(), threshold);
    BOOST_REQUIRE_EQUAL(scaling,               parsedInput.scaling());
    BOOST_REQUIRE_EQUAL(radiiSet,              parsedInput.radiiSet());
    BOOST_REQUIRE_EQUAL(mode,                  parsedInput.mode());
    for (size_t i = 0; i < spheres.size(); ++i) {
	    for (size_t j = 0; j < 3; ++j) {
	    	BOOST_REQUIRE_CLOSE(spheres[i].center(j),      parsedInput.spheres(i).center(j), threshold);
	    }
	    BOOST_REQUIRE_CLOSE(spheres[i].radius(),      parsedInput.spheres(i).radius(), threshold);
    }
    BOOST_REQUIRE_EQUAL(solverType,            parsedInput.solverType());
    BOOST_REQUIRE_EQUAL(equationType,          parsedInput.equationType());
    BOOST_REQUIRE_CLOSE(probeRadius,           parsedInput.probeRadius(), 1.0e-10);
    BOOST_REQUIRE_EQUAL(greenInsideType,       parsedInput.greenInsideType());
    BOOST_REQUIRE_EQUAL(greenOutsideType,      parsedInput.greenOutsideType());
    BOOST_REQUIRE_EQUAL(derivativeInsideType,  parsedInput.derivativeInsideType());
    BOOST_REQUIRE_EQUAL(derivativeOutsideType, parsedInput.derivativeOutsideType());
    BOOST_REQUIRE_CLOSE(epsilonInside,         parsedInput.epsilonInside(), threshold);
    BOOST_REQUIRE_CLOSE(epsilonOutside,        parsedInput.epsilonOutside(), threshold);
}

BOOST_AUTO_TEST_SUITE_END()
