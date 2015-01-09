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

#define BOOST_TEST_MODULE GePolCavity

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <vector>
#include <cmath>

#include "Config.hpp"

#include <Eigen/Dense>

#include <boost/filesystem.hpp>

#include "GePolCavity.hpp"
#include "PhysicalConstants.hpp"
#include "TestingMolecules.hpp"

namespace fs = boost::filesystem;

struct GePolCavityC1Test {
    GePolCavity cavity;
    GePolCavityC1Test() {  SetUp(); }
    void SetUp() {
	Molecule point = dummy<0>();
        double area = 0.4;
        double probeRadius = 0.0;
        double minRadius = 100.0;
        cavity = GePolCavity(point, area, probeRadius, minRadius);
        fs::rename("PEDRA.OUT", "PEDRA.OUT.c1");
        fs::rename("cavity.off", "cavity.off.c1");
    }
};

BOOST_FIXTURE_TEST_SUITE(GePolCavityC1, GePolCavityC1Test)

/*! \class GePolCavity
 *  \test \b GePolCavityC1Test_size tests GePol cavity size for a point charge in C1 symmetry with added spheres
 */
BOOST_FIXTURE_TEST_CASE(size, GePolCavityC1Test)
{
    int size = 32;
    int actualSize = cavity.size();
    BOOST_REQUIRE_EQUAL(size, actualSize);
}

/*! \class GePolCavity
 *  \test \b GePolCavityC1Test_irreducible_size tests GePol cavity irreducible size for a point charge in C1 symmetry with added spheres
 */
BOOST_FIXTURE_TEST_CASE(irreducible_size, GePolCavityC1Test)
{
    int size = 32;
    int actualSize = cavity.irreducible_size();
    BOOST_REQUIRE_EQUAL(size, actualSize);
}

/*! \class GePolCavity
 *  \test \b GePolCavityC1Test_area tests GePol cavity surface area for a point charge in C1 symmetry with added spheres
 */
BOOST_FIXTURE_TEST_CASE(area, GePolCavityC1Test)
{
    double area = 4.0 * M_PI * pow(1.0, 2);
    double actualArea = cavity.elementArea().sum();
    BOOST_REQUIRE_CLOSE(area, actualArea, 1.0e-12);
}

/*! \class GePolCavity
 *  \test \b GePolCavityC1Test_volume tests GePol cavity volume for a point charge in C1 symmetry with added spheres
 */
BOOST_FIXTURE_TEST_CASE(volume, GePolCavityC1Test)
{
    double volume = 4.0 * M_PI * pow(1.0, 3) / 3.0;
    Eigen::Matrix3Xd elementCenter = cavity.elementCenter();
    Eigen::Matrix3Xd elementNormal = cavity.elementNormal();
    double actualVolume = 0;
    for ( int i = 0; i < cavity.size(); ++i ) {
        actualVolume += cavity.elementArea(i) * elementCenter.col(i).dot(elementNormal.col(
                            i));
    }
    actualVolume /= 3;
    BOOST_REQUIRE_CLOSE(volume, actualVolume, 1.0e-12);
}

BOOST_AUTO_TEST_SUITE_END()

struct GePolCavityC2Test {
    GePolCavity cavity;
    GePolCavityC2Test() { SetUp(); }
    void SetUp() {
	Molecule point = dummy<1>();
        double area = 0.4;
        double probeRadius = 0.0;
        double minRadius = 100.0;
        cavity = GePolCavity(point, area, probeRadius, minRadius);
        fs::rename("PEDRA.OUT", "PEDRA.OUT.c2");
        fs::rename("cavity.off", "cavity.off.c2");
    }
};

BOOST_FIXTURE_TEST_SUITE(GePolCavityC2, GePolCavityC2Test)

/*! \class GePolCavity
 *  \test \b GePolCavityC2Test_size tests GePol cavity size for a point charge in C2 symmetry
 */
BOOST_FIXTURE_TEST_CASE(size, GePolCavityC2Test)
{
    int size = 32;
    int actualSize = cavity.size();
    BOOST_REQUIRE_EQUAL(size, actualSize);
}

/*! \class GePolCavity
 *  \test \b GePolCavityC2Test_irreducible_size tests GePol cavity irreducible size for a point charge in C2 symmetry
 */
BOOST_FIXTURE_TEST_CASE(irreducible_size, GePolCavityC2Test)
{
    int size = 16;
    int actualSize = cavity.irreducible_size();
    BOOST_REQUIRE_EQUAL(size, actualSize);
}

/*! \class GePolCavity
 *  \test \b GePolCavityC2Test_area tests GePol cavity surface area for a point charge in C2 symmetry
 */
BOOST_FIXTURE_TEST_CASE(area, GePolCavityC2Test)
{
    double area = 4.0 * M_PI * pow(1.0, 2);
    double actualArea = cavity.elementArea().sum();
    BOOST_REQUIRE_CLOSE(area, actualArea, 1.0e-12);
}

/*! \class GePolCavity
 *  \test \b GePolCavityC2Test_volume tests GePol cavity volume for a point charge in C2 symmetry
 */
BOOST_FIXTURE_TEST_CASE(volume, GePolCavityC2Test)
{
    double volume = 4.0 * M_PI * pow(1.0, 3) / 3.0;
    Eigen::Matrix3Xd elementCenter = cavity.elementCenter();
    Eigen::Matrix3Xd elementNormal = cavity.elementNormal();
    double actualVolume = 0;
    for ( int i = 0; i < cavity.size(); ++i ) {
        actualVolume += cavity.elementArea(i) * elementCenter.col(i).dot(elementNormal.col(
                            i));
    }
    actualVolume /= 3;
    BOOST_REQUIRE_CLOSE(volume, actualVolume, 1.0e-12);
}

BOOST_AUTO_TEST_SUITE_END()

struct GePolCavityCsTest {
    GePolCavity cavity;
    GePolCavityCsTest() { SetUp(); }
    void SetUp() {
	Molecule point = dummy<2>();
        double area = 0.4;
        double probeRadius = 0.0;
        double minRadius = 100.0;
        cavity = GePolCavity(point, area, probeRadius, minRadius);
        fs::rename("PEDRA.OUT", "PEDRA.OUT.cs");
        fs::rename("cavity.off", "cavity.off.cs");
    }
};

BOOST_FIXTURE_TEST_SUITE(GePolCavityCs, GePolCavityCsTest)

/*! \class GePolCavity
 *  \test \b GePolCavityCsTest_size tests GePol cavity size for a point charge in Cs symmetry with added spheres
 */
BOOST_FIXTURE_TEST_CASE(size, GePolCavityCsTest)
{
    int size = 32;
    int actualSize = cavity.size();
    BOOST_REQUIRE_EQUAL(size, actualSize);
}

/*! \class GePolCavity
 *  \test \b GePolCavityCsTest_irreducible_size tests GePol cavity irreducible size for a point charge in Cs symmetry
 */
BOOST_FIXTURE_TEST_CASE(irreducible_size, GePolCavityCsTest)
{
    int size = 16;
    int actualSize = cavity.irreducible_size();
    BOOST_REQUIRE_EQUAL(size, actualSize);
}

/*! \class GePolCavity
 *  \test \b GePolCavityCsTest_area tests GePol cavity surface area for a point charge in Cs symmetry
 */
BOOST_FIXTURE_TEST_CASE(area, GePolCavityCsTest)
{
    double area = 4.0 * M_PI * pow(1.0, 2);
    double actualArea = cavity.elementArea().sum();
    BOOST_REQUIRE_CLOSE(area, actualArea, 1.0e-12);
}

/*! \class GePolCavity
 *  \test \b GePolCavityCsTest_volume tests GePol cavity volume for a point charge in Cs symmetry
 */
BOOST_FIXTURE_TEST_CASE(volume, GePolCavityCsTest)
{
    double volume = 4.0 * M_PI * pow(1.0, 3) / 3.0;
    Eigen::Matrix3Xd elementCenter = cavity.elementCenter();
    Eigen::Matrix3Xd elementNormal = cavity.elementNormal();
    double actualVolume = 0;
    for ( int i = 0; i < cavity.size(); ++i ) {
        actualVolume += cavity.elementArea(i) * elementCenter.col(i).dot(elementNormal.col(
                            i));
    }
    actualVolume /= 3;
    BOOST_REQUIRE_CLOSE(volume, actualVolume, 1.0e-12);
}

BOOST_AUTO_TEST_SUITE_END()

struct GePolCavityCiTest {
    GePolCavity cavity;
    GePolCavityCiTest() { SetUp(); }
    void SetUp() {
	Molecule point = dummy<3>();
        double area = 0.4;
        double probeRadius = 0.0;
        double minRadius = 100.0;
        cavity = GePolCavity(point, area, probeRadius, minRadius);
        fs::rename("PEDRA.OUT", "PEDRA.OUT.ci");
        fs::rename("cavity.off", "cavity.off.ci");
    }
};

BOOST_FIXTURE_TEST_SUITE(GePolCavityCi, GePolCavityCiTest)

/*! \class GePolCavity
 *  \test \b GePolCavityCiTest_size tests GePol cavity size for a point charge in Ci symmetry
 */
BOOST_FIXTURE_TEST_CASE(size, GePolCavityCiTest)
{
    int size = 32;
    int actualSize = cavity.size();
    BOOST_REQUIRE_EQUAL(size, actualSize);
}

/*! \class GePolCavity
 *  \test \b GePolCavityCiTest_irreducible_size tests GePol cavity irreducible size for a point charge in Ci symmetry
 */
BOOST_FIXTURE_TEST_CASE(irreducible_size, GePolCavityCiTest)
{
    int size = 16;
    int actualSize = cavity.irreducible_size();
    BOOST_REQUIRE_EQUAL(size, actualSize);
}

/*! \class GePolCavity
 *  \test \b GePolCavityCiTest_area tests GePol cavity surface area for a point charge in Ci symmetry
 */
BOOST_FIXTURE_TEST_CASE(area, GePolCavityCiTest)
{
    double area = 4.0 * M_PI * pow(1.0, 2);
    double actualArea = cavity.elementArea().sum();
    BOOST_REQUIRE_CLOSE(area, actualArea, 1.0e-12);
}

/*! \class GePolCavity
 *  \test \b GePolCavityCiTest_volume tests GePol cavity volume for a point charge in Ci symmetry
 */
BOOST_FIXTURE_TEST_CASE(volume, GePolCavityCiTest)
{
    double volume = 4.0 * M_PI * pow(1.0, 3) / 3.0;
    Eigen::Matrix3Xd elementCenter = cavity.elementCenter();
    Eigen::Matrix3Xd elementNormal = cavity.elementNormal();
    double actualVolume = 0;
    for ( int i = 0; i < cavity.size(); ++i ) {
        actualVolume += cavity.elementArea(i) * elementCenter.col(i).dot(elementNormal.col(
                            i));
    }
    actualVolume /= 3;
    BOOST_REQUIRE_CLOSE(volume, actualVolume, 1.0e-12);
}

BOOST_AUTO_TEST_SUITE_END()

struct GePolCavityD2Test {
    GePolCavity cavity;
    GePolCavityD2Test() { SetUp(); }
    void SetUp() {
	Molecule point = dummy<4>();
        double area = 0.4;
        double probeRadius = 0.0;
        double minRadius = 100.0;
        cavity = GePolCavity(point, area, probeRadius, minRadius);
        fs::rename("PEDRA.OUT", "PEDRA.OUT.d2");
        fs::rename("cavity.off", "cavity.off.d2");
    }
};

BOOST_FIXTURE_TEST_SUITE(GePolCavityD2, GePolCavityD2Test)

/*! \class GePolCavity
 *  \test \b GePolCavityD2Test_size tests GePol cavity size for a point charge in D2 symmetry
 */
BOOST_FIXTURE_TEST_CASE(size, GePolCavityD2Test)
{
    int size = 32;
    int actualSize = cavity.size();
    BOOST_REQUIRE_EQUAL(size, actualSize);
}

/*! \class GePolCavity
 *  \test \b GePolCavityD2Test_irreducible_size tests GePol cavity irreducible size for a point charge in D2 symmetry
 */
BOOST_FIXTURE_TEST_CASE(irreducible_size, GePolCavityD2Test)
{
    int size = 8;
    int actualSize = cavity.irreducible_size();
    BOOST_REQUIRE_EQUAL(size, actualSize);
}

/*! \class GePolCavity
 *  \test \b GePolCavityD2Test_area tests GePol cavity surface area for a point charge in D2 symmetry
 */
BOOST_FIXTURE_TEST_CASE(area, GePolCavityD2Test)
{
    double area = 4.0 * M_PI * pow(1.0, 2);
    double actualArea = cavity.elementArea().sum();
    BOOST_REQUIRE_CLOSE(area, actualArea, 1.0e-12);
}

/*! \class GePolCavity
 *  \test \b GePolCavityD2Test_volume tests GePol cavity volume for a point charge in D2 symmetry
 */
BOOST_FIXTURE_TEST_CASE(volume, GePolCavityD2Test)
{
    double volume = 4.0 * M_PI * pow(1.0, 3) / 3.0;
    Eigen::Matrix3Xd elementCenter = cavity.elementCenter();
    Eigen::Matrix3Xd elementNormal = cavity.elementNormal();
    double actualVolume = 0;
    for ( int i = 0; i < cavity.size(); ++i ) {
        actualVolume += cavity.elementArea(i) * elementCenter.col(i).dot(elementNormal.col(
                            i));
    }
    actualVolume /= 3;
    BOOST_REQUIRE_CLOSE(volume, actualVolume, 1.0e-12);
}

BOOST_AUTO_TEST_SUITE_END()

struct GePolCavityC2vTest {
    GePolCavity cavity;
    GePolCavityC2vTest() { SetUp(); }
    void SetUp() {
	Molecule point = dummy<5>();
        double area = 0.4;
        double probeRadius = 0.0;
        double minRadius = 100.0;
        cavity = GePolCavity(point, area, probeRadius, minRadius);
        fs::rename("PEDRA.OUT", "PEDRA.OUT.c2v");
        fs::rename("cavity.off", "cavity.off.c2v");
    }
};

BOOST_FIXTURE_TEST_SUITE(GePolCavityC2v, GePolCavityC2vTest)

/*! \class GePolCavity
 *  \test \b GePolCavityC2vTest_size tests GePol cavity size for a point charge in C2v symmetry
 */
BOOST_FIXTURE_TEST_CASE(size, GePolCavityC2vTest)
{
    int size = 32;
    int actualSize = cavity.size();
    BOOST_REQUIRE_EQUAL(size, actualSize);
}

/*! \class GePolCavity
 *  \test \b GePolCavityC2vTest_irreducible_size tests GePol cavity irreducible size for a point charge in C2v symmetry
 */
BOOST_FIXTURE_TEST_CASE(irreducible_size, GePolCavityC2vTest)
{
    int size = 8;
    int actualSize = cavity.irreducible_size();
    BOOST_REQUIRE_EQUAL(size, actualSize);
}

/*! \class GePolCavity
 *  \test \b GePolCavityC2vTest_area tests GePol cavity surface area for a point charge in C2v symmetry
 */
BOOST_FIXTURE_TEST_CASE(area, GePolCavityC2vTest)
{
    double area = 4.0 * M_PI * pow(1.0, 2);
    double actualArea = cavity.elementArea().sum();
    BOOST_REQUIRE_CLOSE(area, actualArea, 1.0e-12);
}

/*! \class GePolCavity
 *  \test \b GePolCavityC2vTest_volume tests GePol cavity volume for a point charge in C2v symmetry
 */
BOOST_FIXTURE_TEST_CASE(volume, GePolCavityC2vTest)
{
    double volume = 4.0 * M_PI * pow(1.0, 3) / 3.0;
    Eigen::Matrix3Xd elementCenter = cavity.elementCenter();
    Eigen::Matrix3Xd elementNormal = cavity.elementNormal();
    double actualVolume = 0;
    for ( int i = 0; i < cavity.size(); ++i ) {
        actualVolume += cavity.elementArea(i) * elementCenter.col(i).dot(elementNormal.col(
                            i));
    }
    actualVolume /= 3;
    BOOST_REQUIRE_CLOSE(volume, actualVolume, 1.0e-12);
}

BOOST_AUTO_TEST_SUITE_END()

struct GePolCavityC2hTest {
    GePolCavity cavity;
    GePolCavityC2hTest() { SetUp(); }
    void SetUp() {
	Molecule point = dummy<6>();
        double area = 0.4;
        double probeRadius = 0.0;
        double minRadius = 100.0;
        cavity = GePolCavity(point, area, probeRadius, minRadius);
        fs::rename("PEDRA.OUT", "PEDRA.OUT.c2h");
        fs::rename("cavity.off", "cavity.off.c2h");
    }
};

BOOST_FIXTURE_TEST_SUITE(GePolCavityC2h, GePolCavityC2hTest)

/*! \class GePolCavity
 *  \test \b GePolCavityC2hTest_size tests GePol cavity size for a point charge in C2h symmetry
 */
BOOST_FIXTURE_TEST_CASE(size, GePolCavityC2hTest)
{
    int size = 32;
    int actualSize = cavity.size();
    BOOST_REQUIRE_EQUAL(size, actualSize);
}

/*! \class GePolCavity
 *  \test \b GePolCavityC2hTest_irreducible_size tests GePol cavity irreducible size for a point charge in C2h symmetry
 */
BOOST_FIXTURE_TEST_CASE(irreducible_size, GePolCavityC2hTest)
{
    int size = 8;
    int actualSize = cavity.irreducible_size();
    BOOST_REQUIRE_EQUAL(size, actualSize);
}

/*! \class GePolCavity
 *  \test \b GePolCavityC2hTest_area tests GePol cavity surface area for a point charge in C2h symmetry
 */
BOOST_FIXTURE_TEST_CASE(area, GePolCavityC2hTest)
{
    double area = 4.0 * M_PI * pow(1.0, 2);
    double actualArea = cavity.elementArea().sum();
    BOOST_REQUIRE_CLOSE(area, actualArea, 1.0e-12);
}

/*! \class GePolCavity
 *  \test \b GePolCavityC2hTest_volume tests GePol cavity volume for a point charge in C2h symmetry
 */
BOOST_FIXTURE_TEST_CASE(volume, GePolCavityC2hTest)
{
    double volume = 4.0 * M_PI * pow(1.0, 3) / 3.0;
    Eigen::Matrix3Xd elementCenter = cavity.elementCenter();
    Eigen::Matrix3Xd elementNormal = cavity.elementNormal();
    double actualVolume = 0;
    for ( int i = 0; i < cavity.size(); ++i ) {
        actualVolume += cavity.elementArea(i) * elementCenter.col(i).dot(elementNormal.col(
                            i));
    }
    actualVolume /= 3;
    BOOST_REQUIRE_CLOSE(volume, actualVolume, 1.0e-12);
}

BOOST_AUTO_TEST_SUITE_END()

struct GePolCavityD2hTest {
    GePolCavity cavity;
    GePolCavityD2hTest() { SetUp(); }
    void SetUp() {
	Molecule point = dummy<7>();
        double area = 0.4;
        double probeRadius = 0.0;
        double minRadius = 100.0;
        cavity = GePolCavity(point, area, probeRadius, minRadius);
        fs::rename("PEDRA.OUT", "PEDRA.OUT.d2h");
        fs::rename("cavity.off", "cavity.off.d2h");
    }
};

BOOST_FIXTURE_TEST_SUITE(GePolCavityD2h, GePolCavityD2hTest)

/*! \class GePolCavity
 *  \test \b GePolCavityD2hTest_size tests GePol cavity size for a point charge in D2h symmetry
 */
BOOST_FIXTURE_TEST_CASE(size, GePolCavityD2hTest)
{
    int size = 32;
    int actualSize = cavity.size();
    BOOST_REQUIRE_EQUAL(size, actualSize);
}

/*! \class GePolCavity
 *  \test \b GePolCavityD2hTest_irreducible_size tests GePol cavity irreducible size for a point charge in D2h symmetry
 */
BOOST_FIXTURE_TEST_CASE(irreducible_size, GePolCavityD2hTest)
{
    int size = 4;
    int actualSize = cavity.irreducible_size();
    BOOST_REQUIRE_EQUAL(size, actualSize);
}

/*! \class GePolCavity
 *  \test \b GePolCavityD2hTest_area tests GePol cavity surface area for a point charge in D2h symmetry
 */
BOOST_FIXTURE_TEST_CASE(area, GePolCavityD2hTest)
{
    double area = 4.0 * M_PI * pow(1.0, 2);
    double actualArea = cavity.elementArea().sum();
    BOOST_REQUIRE_CLOSE(area, actualArea, 1.0e-12);
}

/*! \class GePolCavity
 *  \test \b GePolCavityD2hTest_volume tests GePol cavity volume for a point charge in D2h symmetry
 */
BOOST_FIXTURE_TEST_CASE(volume, GePolCavityD2hTest)
{
    double volume = 4.0 * M_PI * pow(1.0, 3) / 3.0;
    Eigen::Matrix3Xd elementCenter = cavity.elementCenter();
    Eigen::Matrix3Xd elementNormal = cavity.elementNormal();
    double actualVolume = 0;
    for ( int i = 0; i < cavity.size(); ++i ) {
        actualVolume += cavity.elementArea(i) * elementCenter.col(i).dot(elementNormal.col(
                            i));
    }
    actualVolume /= 3;
    BOOST_REQUIRE_CLOSE(volume, actualVolume, 1.0e-12);
}

BOOST_AUTO_TEST_SUITE_END()
