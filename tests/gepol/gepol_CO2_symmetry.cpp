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

#define BOOST_TEST_MODULE GePolCavityCO2

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <vector>
#include <cmath>

#include "Config.hpp"

#include <Eigen/Dense>

#include <boost/filesystem.hpp>

#include "GePolCavity.hpp"
#include "PhysicalConstants.hpp"
#include "Symmetry.hpp"

namespace fs = boost::filesystem;

struct GePolCavityCO2C1Test {
    GePolCavity cavity;
    GePolCavityCO2C1Test() { SetUp(); }
    void SetUp() {
        Eigen::Vector3d	C1( 0.0000000000, 0.0000000000, 0.0000000000);
        Eigen::Vector3d	O1( 2.1316110791, 0.0000000000, 0.0000000000);
        Eigen::Vector3d	O2(-2.1316110791, 0.0000000000, 0.0000000000);
        std::vector<Sphere> spheres;
        double radiusC = (1.70 * 1.20) / convertBohrToAngstrom;
        double radiusO = (1.52 * 1.20) / convertBohrToAngstrom;
        Sphere sph1(C1, radiusC);
        Sphere sph2(O1, radiusO);
        Sphere sph3(O2, radiusO);
        spheres.push_back(sph1);
        spheres.push_back(sph2);
        spheres.push_back(sph3);
        double area = 0.2 / convertBohr2ToAngstrom2;
        double probeRadius = 1.385 / convertBohrToAngstrom;
        double minRadius = 0.2 / convertBohrToAngstrom;
        // C1
        Symmetry pGroup = buildGroup(0, 0, 0, 0);
        cavity = GePolCavity(spheres, area, probeRadius, minRadius, pGroup);
        fs::rename("PEDRA.OUT", "PEDRA.co2.c1");
        fs::rename("cavity.off", "cavity.co2.c1");
    }
};

BOOST_FIXTURE_TEST_SUITE(GePolCavityCO2C1, GePolCavityCO2C1Test)

/*! \class GePolCavity
 *  \test \b GePolCavityCO2C1Test_size tests GePol cavity size for CO2 in C1 symmetry
 */
BOOST_FIXTURE_TEST_CASE(size, GePolCavityCO2C1Test)
{
    int size = 448;
    int actualSize = cavity.size();
    BOOST_REQUIRE_EQUAL(size, actualSize);
}

/*! \class GePolCavity
 *  \test \b GePolCavityCO2C1Test_irreducible_size tests GePol cavity irreducible size for CO2 in C1 symmetry
 */
BOOST_FIXTURE_TEST_CASE(irreducible_size, GePolCavityCO2C1Test)
{
    int size = 448;
    int actualSize = cavity.irreducible_size();
    BOOST_REQUIRE_EQUAL(size, actualSize);
}

/*! \class GePolCavity
 *  \test \b GePolCavityCO2C1Test_area tests GePol cavity surface area for CO2 in C1 symmetry
 */
BOOST_FIXTURE_TEST_CASE(area, GePolCavityCO2C1Test)
{
    double area = 250.68176442433020;
    double actualArea = cavity.elementArea().sum();
    BOOST_REQUIRE_CLOSE(area, actualArea, 1.0e-10);
}

/*! \class GePolCavity
 *  \test \b GePolCavityCO2C1Test_volume tests GePol cavity volume for CO2 in C1 symmetry
 */
BOOST_FIXTURE_TEST_CASE(volume, GePolCavityCO2C1Test)
{
    double volume = 352.55869984340751;
    Eigen::Matrix3Xd elementCenter = cavity.elementCenter();
    Eigen::Matrix3Xd elementNormal = cavity.elementNormal();
    double actualVolume = 0;
    for ( int i = 0; i < cavity.size(); ++i ) {
        actualVolume += cavity.elementArea(i) * elementCenter.col(i).dot(elementNormal.col(
                            i));
    }
    actualVolume /= 3;
    BOOST_REQUIRE_CLOSE(volume, actualVolume, 1.0e-10);
}

BOOST_AUTO_TEST_SUITE_END()

struct GePolCavityCO2CsTest {
    GePolCavity cavity;
    GePolCavityCO2CsTest() { SetUp(); }
    void SetUp() {
        Eigen::Vector3d	C1( 0.0000000000, 0.0000000000, 0.0000000000);
        Eigen::Vector3d	O1( 2.1316110791, 0.0000000000, 0.0000000000);
        Eigen::Vector3d	O2(-2.1316110791, 0.0000000000, 0.0000000000);
        std::vector<Sphere> spheres;
        double radiusC = (1.70 * 1.20) / convertBohrToAngstrom;
        double radiusO = (1.52 * 1.20) / convertBohrToAngstrom;
        Sphere sph1(C1, radiusC);
        Sphere sph2(O1, radiusO);
        Sphere sph3(O2, radiusO);
        spheres.push_back(sph1);
        spheres.push_back(sph2);
        spheres.push_back(sph3);
        double area = 0.2 / convertBohr2ToAngstrom2;
        double probeRadius = 1.385 / convertBohrToAngstrom;
        double minRadius = 0.2 / convertBohrToAngstrom;
        // Cs as generated by Oyz
        Symmetry pGroup = buildGroup(1, 1, 0, 0);
        cavity = GePolCavity(spheres, area, probeRadius, minRadius, pGroup);
        fs::rename("PEDRA.OUT", "PEDRA.co2.cs");
        fs::rename("cavity.off", "cavity.co2.cs");
    }
};

BOOST_FIXTURE_TEST_SUITE(GePolCavityCO2Cs, GePolCavityCO2CsTest)

/*! \class GePolCavity
 *  \test \b GePolCavityCO2CsTest_size tests GePol cavity size for CO2 in Cs symmetry
 */
BOOST_FIXTURE_TEST_CASE(size, GePolCavityCO2CsTest)
{
    int size = 448;
    int actualSize = cavity.size();
    BOOST_REQUIRE_EQUAL(size, actualSize);
}

/*! \class GePolCavity
 *  \test \b GePolCavityCO2CsTest_irreducible_size tests GePol cavity irreducible size for CO2 in Cs symmetry
 */
BOOST_FIXTURE_TEST_CASE(irreducible_size, GePolCavityCO2CsTest)
{
    int size = 224;
    int actualSize = cavity.irreducible_size();
    BOOST_REQUIRE_EQUAL(size, actualSize);
}

/*! \class GePolCavity
 *  \test \b GePolCavityCO2CsTest_area tests GePol cavity surface area for CO2 in Cs symmetry
 */
BOOST_FIXTURE_TEST_CASE(area, GePolCavityCO2CsTest)
{
    double area = 250.68176442433020;
    double actualArea = cavity.elementArea().sum();
    BOOST_REQUIRE_CLOSE(area, actualArea, 1.0e-10);
}

/*! \class GePolCavity
 *  \test \b GePolCavityCO2CsTest_volume tests GePol cavity volume for CO2 in Cs symmetry
 */
BOOST_FIXTURE_TEST_CASE(volume, GePolCavityCO2CsTest)
{
    double volume = 352.55869984340751;
    Eigen::Matrix3Xd elementCenter = cavity.elementCenter();
    Eigen::Matrix3Xd elementNormal = cavity.elementNormal();
    double actualVolume = 0;
    for ( int i = 0; i < cavity.size(); ++i ) {
        actualVolume += cavity.elementArea(i) * elementCenter.col(i).dot(elementNormal.col(
                            i));
    }
    actualVolume /= 3;
    BOOST_REQUIRE_CLOSE(volume, actualVolume, 1.0e-10);
}

BOOST_AUTO_TEST_SUITE_END()

struct GePolCavityCO2C2Test {
    GePolCavity cavity;
    GePolCavityCO2C2Test() { SetUp(); }
    void SetUp() {
        Eigen::Vector3d	C1( 0.0000000000, 0.0000000000, 0.0000000000);
        Eigen::Vector3d	O1( 2.1316110791, 0.0000000000, 0.0000000000);
        Eigen::Vector3d	O2(-2.1316110791, 0.0000000000, 0.0000000000);
        std::vector<Sphere> spheres;
        double radiusC = (1.70 * 1.20) / convertBohrToAngstrom;
        double radiusO = (1.52 * 1.20) / convertBohrToAngstrom;
        Sphere sph1(C1, radiusC);
        Sphere sph2(O1, radiusO);
        Sphere sph3(O2, radiusO);
        spheres.push_back(sph1);
        spheres.push_back(sph2);
        spheres.push_back(sph3);
        double area = 0.2 / convertBohr2ToAngstrom2;
        double probeRadius = 1.385 / convertBohrToAngstrom;
        double minRadius = 0.2 / convertBohrToAngstrom;
        // C2 as generated by C2z
        Symmetry pGroup = buildGroup(1, 3, 0, 0);
        cavity = GePolCavity(spheres, area, probeRadius, minRadius, pGroup);
        fs::rename("PEDRA.OUT", "PEDRA.co2.c2");
        fs::rename("cavity.off", "cavity.co2.c2");
    }
};

BOOST_FIXTURE_TEST_SUITE(GePolCavityCO2C2, GePolCavityCO2C2Test)

/*! \class GePolCavity
 *  \test \b GePolCavityCO2C2Test_size tests GePol cavity size for CO2 in C2 symmetry
 */
BOOST_FIXTURE_TEST_CASE(size, GePolCavityCO2C2Test)
{
    int size = 448;
    int actualSize = cavity.size();
    BOOST_REQUIRE_EQUAL(size, actualSize);
}

/*! \class GePolCavity
 *  \test \b GePolCavityCO2C2Test_irreducible_size tests GePol cavity irreducible size for CO2 in C2 symmetry
 */
BOOST_FIXTURE_TEST_CASE(irreducible_size, GePolCavityCO2C2Test)
{
    int size = 224;
    int actualSize = cavity.irreducible_size();
    BOOST_REQUIRE_EQUAL(size, actualSize);
}

/*! \class GePolCavity
 *  \test \b GePolCavityCO2C2Test_area tests GePol cavity surface area for CO2 in C2 symmetry
 */
BOOST_FIXTURE_TEST_CASE(area, GePolCavityCO2C2Test)
{
    double area = 250.68176442433020;
    double actualArea = cavity.elementArea().sum();
    BOOST_REQUIRE_CLOSE(area, actualArea, 1.0e-10);
}

/*! \class GePolCavity
 *  \test \b GePolCavityCO2C2Test_volume tests GePol cavity volume for CO2 in C2 symmetry
 */
BOOST_FIXTURE_TEST_CASE(volume, GePolCavityCO2C2Test)
{
    double volume = 352.55869984340751;
    Eigen::Matrix3Xd elementCenter = cavity.elementCenter();
    Eigen::Matrix3Xd elementNormal = cavity.elementNormal();
    double actualVolume = 0;
    for ( int i = 0; i < cavity.size(); ++i ) {
        actualVolume += cavity.elementArea(i) * elementCenter.col(i).dot(elementNormal.col(
                            i));
    }
    actualVolume /= 3;
    BOOST_REQUIRE_CLOSE(volume, actualVolume, 1.0e-10);
}

BOOST_AUTO_TEST_SUITE_END()

struct GePolCavityCO2CiTest {
    GePolCavity cavity;
    GePolCavityCO2CiTest() { SetUp(); }
    void SetUp() {
        Eigen::Vector3d	C1( 0.0000000000, 0.0000000000, 0.0000000000);
        Eigen::Vector3d	O1( 2.1316110791, 0.0000000000, 0.0000000000);
        Eigen::Vector3d	O2(-2.1316110791, 0.0000000000, 0.0000000000);
        std::vector<Sphere> spheres;
        double radiusC = (1.70 * 1.20) / convertBohrToAngstrom;
        double radiusO = (1.52 * 1.20) / convertBohrToAngstrom;
        Sphere sph1(C1, radiusC);
        Sphere sph2(O1, radiusO);
        Sphere sph3(O2, radiusO);
        spheres.push_back(sph1);
        spheres.push_back(sph2);
        spheres.push_back(sph3);
        double area = 0.2 / convertBohr2ToAngstrom2;
        double probeRadius = 1.385 / convertBohrToAngstrom;
        double minRadius = 0.2 / convertBohrToAngstrom;
        // Ci as generated by i
        Symmetry pGroup = buildGroup(1, 7, 0, 0);
        cavity = GePolCavity(spheres, area, probeRadius, minRadius, pGroup);
        fs::rename("PEDRA.OUT", "PEDRA.co2.ci");
        fs::rename("cavity.off", "cavity.co2.ci");
    }
};

BOOST_FIXTURE_TEST_SUITE(GePolCavityCO2Ci, GePolCavityCO2CiTest)

/*! \class GePolCavity
 *  \test \b GePolCavityCO2CiTest_size tests GePol cavity size for CO2 in Ci symmetry
 */
BOOST_FIXTURE_TEST_CASE(size, GePolCavityCO2CiTest)
{
    int size = 448;
    int actualSize = cavity.size();
    BOOST_REQUIRE_EQUAL(size, actualSize);
}

/*! \class GePolCavity
 *  \test \b GePolCavityCO2CiTest_irreducible_size tests GePol cavity irreducible size for CO2 in Ci symmetry
 */
BOOST_FIXTURE_TEST_CASE(irreducible_size, GePolCavityCO2CiTest)
{
    int size = 224;
    int actualSize = cavity.irreducible_size();
    BOOST_REQUIRE_EQUAL(size, actualSize);
}

/*! \class GePolCavity
 *  \test \b GePolCavityCO2CiTest_area tests GePol cavity surface area for CO2 in Ci symmetry
 */
BOOST_FIXTURE_TEST_CASE(area, GePolCavityCO2CiTest)
{
    double area = 250.68176442433020;
    double actualArea = cavity.elementArea().sum();
    BOOST_REQUIRE_CLOSE(area, actualArea, 1.0e-10);
}

/*! \class GePolCavity
 *  \test \b GePolCavityCO2CiTest_volume tests GePol cavity volume for CO2 in Ci symmetry
 */
BOOST_FIXTURE_TEST_CASE(volume, GePolCavityCO2CiTest)
{
    double volume = 352.55869984340751;
    Eigen::Matrix3Xd elementCenter = cavity.elementCenter();
    Eigen::Matrix3Xd elementNormal = cavity.elementNormal();
    double actualVolume = 0;
    for ( int i = 0; i < cavity.size(); ++i ) {
        actualVolume += cavity.elementArea(i) * elementCenter.col(i).dot(elementNormal.col(
                            i));
    }
    actualVolume /= 3;
    BOOST_REQUIRE_CLOSE(volume, actualVolume, 1.0e-10);
}

BOOST_AUTO_TEST_SUITE_END()

struct GePolCavityCO2C2hTest {
    GePolCavity cavity;
    GePolCavityCO2C2hTest() { SetUp(); }
    void SetUp() {
        Eigen::Vector3d	C1( 0.0000000000, 0.0000000000, 0.0000000000);
        Eigen::Vector3d	O1( 2.1316110791, 0.0000000000, 0.0000000000);
        Eigen::Vector3d	O2(-2.1316110791, 0.0000000000, 0.0000000000);
        std::vector<Sphere> spheres;
        double radiusC = (1.70 * 1.20) / convertBohrToAngstrom;
        double radiusO = (1.52 * 1.20) / convertBohrToAngstrom;
        Sphere sph1(C1, radiusC);
        Sphere sph2(O1, radiusO);
        Sphere sph3(O2, radiusO);
        spheres.push_back(sph1);
        spheres.push_back(sph2);
        spheres.push_back(sph3);
        double area = 0.2 / convertBohr2ToAngstrom2;
        double probeRadius = 1.385 / convertBohrToAngstrom;
        double minRadius = 0.2 / convertBohrToAngstrom;
        // C2h as generated by Oxy and i
        Symmetry pGroup = buildGroup(2, 4, 7, 0);
        cavity = GePolCavity(spheres, area, probeRadius, minRadius, pGroup);
        fs::rename("PEDRA.OUT", "PEDRA.co2.c2h");
        fs::rename("cavity.off", "cavity.co2.c2h");
    }
};

BOOST_FIXTURE_TEST_SUITE(GePolCavityCO2C2h, GePolCavityCO2C2hTest)

/*! \class GePolCavity
 *  \test \b GePolCavityCO2C2hTest_size tests GePol cavity size for CO2 in C2h symmetry
 */
BOOST_FIXTURE_TEST_CASE(size, GePolCavityCO2C2hTest)
{
    int size = 448;
    int actualSize = cavity.size();
    BOOST_REQUIRE_EQUAL(size, actualSize);
}

/*! \class GePolCavity
 *  \test \b GePolCavityCO2C2hTest_irreducible_size tests GePol cavity irreducible size for CO2 in C2h symmetry
 */
BOOST_FIXTURE_TEST_CASE(irreducible_size, GePolCavityCO2C2hTest)
{
    int size = 112;
    int actualSize = cavity.irreducible_size();
    BOOST_REQUIRE_EQUAL(size, actualSize);
}

/*! \class GePolCavity
 *  \test \b GePolCavityCO2C2hTest_area tests GePol cavity surface area for CO2 in C2h symmetry
 */
BOOST_FIXTURE_TEST_CASE(area, GePolCavityCO2C2hTest)
{
    double area = 250.68176442433020;
    double actualArea = cavity.elementArea().sum();
    BOOST_REQUIRE_CLOSE(area, actualArea, 1.0e-10);
}

/*! \class GePolCavity
 *  \test \b GePolCavityCO2C2hTest_volume tests GePol cavity volume for CO2 in C2h symmetry
 */
BOOST_FIXTURE_TEST_CASE(volume, GePolCavityCO2C2hTest)
{
    double volume = 352.55869984340751;
    Eigen::Matrix3Xd elementCenter = cavity.elementCenter();
    Eigen::Matrix3Xd elementNormal = cavity.elementNormal();
    double actualVolume = 0;
    for ( int i = 0; i < cavity.size(); ++i ) {
        actualVolume += cavity.elementArea(i) * elementCenter.col(i).dot(elementNormal.col(
                            i));
    }
    actualVolume /= 3;
    BOOST_REQUIRE_CLOSE(volume, actualVolume, 1.0e-10);
}

BOOST_AUTO_TEST_SUITE_END()

struct GePolCavityCO2D2Test {
    GePolCavity cavity;
    GePolCavityCO2D2Test() { SetUp(); }
    void SetUp() {
        Eigen::Vector3d	C1( 0.0000000000, 0.0000000000, 0.0000000000);
        Eigen::Vector3d	O1( 2.1316110791, 0.0000000000, 0.0000000000);
        Eigen::Vector3d	O2(-2.1316110791, 0.0000000000, 0.0000000000);
        std::vector<Sphere> spheres;
        double radiusC = (1.70 * 1.20) / convertBohrToAngstrom;
        double radiusO = (1.52 * 1.20) / convertBohrToAngstrom;
        Sphere sph1(C1, radiusC);
        Sphere sph2(O1, radiusO);
        Sphere sph3(O2, radiusO);
        spheres.push_back(sph1);
        spheres.push_back(sph2);
        spheres.push_back(sph3);
        double area = 0.2 / convertBohr2ToAngstrom2;
        double probeRadius = 1.385 / convertBohrToAngstrom;
        double minRadius = 0.2 / convertBohrToAngstrom;
        // D2 as generated by C2z and C2x
        Symmetry pGroup = buildGroup(2, 3, 6, 0);
        cavity = GePolCavity(spheres, area, probeRadius, minRadius, pGroup);
        fs::rename("PEDRA.OUT", "PEDRA.co2.d2");
        fs::rename("cavity.off", "cavity.co2.d2");
    }
};

BOOST_FIXTURE_TEST_SUITE(GePolCavityCO2D2, GePolCavityCO2D2Test)

/*! \class GePolCavity
 *  \test \b GePolCavityCO2D2Test_size tests GePol cavity size for CO2 in D2 symmetry
 */
BOOST_FIXTURE_TEST_CASE(size, GePolCavityCO2D2Test)
{
    int size = 448;
    int actualSize = cavity.size();
    BOOST_REQUIRE_EQUAL(size, actualSize);
}

/*! \class GePolCavity
 *  \test \b GePolCavityCO2D2Test_irreducible_size tests GePol cavity irreducible size for CO2 in D2 symmetry
 */
BOOST_FIXTURE_TEST_CASE(irreducible_size, GePolCavityCO2D2Test)
{
    int size = 112;
    int actualSize = cavity.irreducible_size();
    BOOST_REQUIRE_EQUAL(size, actualSize);
}

/*! \class GePolCavity
 *  \test \b GePolCavityCO2D2Test_area tests GePol cavity surface area for CO2 in D2 symmetry
 */
BOOST_FIXTURE_TEST_CASE(area, GePolCavityCO2D2Test)
{
    double area = 250.68176442433020;
    double actualArea = cavity.elementArea().sum();
    BOOST_REQUIRE_CLOSE(area, actualArea, 1.0e-10);
}

/*! \class GePolCavity
 *  \test \b GePolCavityCO2D2Test_volume tests GePol cavity volume for CO2 in D2 symmetry
 */
BOOST_FIXTURE_TEST_CASE(volume, GePolCavityCO2D2Test)
{
    double volume = 352.55869984340751;
    Eigen::Matrix3Xd elementCenter = cavity.elementCenter();
    Eigen::Matrix3Xd elementNormal = cavity.elementNormal();
    double actualVolume = 0;
    for ( int i = 0; i < cavity.size(); ++i ) {
        actualVolume += cavity.elementArea(i) * elementCenter.col(i).dot(elementNormal.col(
                            i));
    }
    actualVolume /= 3;
    BOOST_REQUIRE_CLOSE(volume, actualVolume, 1.0e-10);
}

BOOST_AUTO_TEST_SUITE_END()

struct GePolCavityCO2C2vTest {
    GePolCavity cavity;
    GePolCavityCO2C2vTest() { SetUp(); }
    void SetUp() {
        Eigen::Vector3d	C1( 0.0000000000, 0.0000000000, 0.0000000000);
        Eigen::Vector3d	O1( 2.1316110791, 0.0000000000, 0.0000000000);
        Eigen::Vector3d	O2(-2.1316110791, 0.0000000000, 0.0000000000);
        std::vector<Sphere> spheres;
        double radiusC = (1.70 * 1.20) / convertBohrToAngstrom;
        double radiusO = (1.52 * 1.20) / convertBohrToAngstrom;
        Sphere sph1(C1, radiusC);
        Sphere sph2(O1, radiusO);
        Sphere sph3(O2, radiusO);
        spheres.push_back(sph1);
        spheres.push_back(sph2);
        spheres.push_back(sph3);
        double area = 0.2 / convertBohr2ToAngstrom2;
        double probeRadius = 1.385 / convertBohrToAngstrom;
        double minRadius = 0.2 / convertBohrToAngstrom;
        // C2v as generated by Oyz and Oxz
        Symmetry pGroup = buildGroup(2, 1, 2, 0);
        cavity = GePolCavity(spheres, area, probeRadius, minRadius, pGroup);
        fs::rename("PEDRA.OUT", "PEDRA.co2.c2v");
        fs::rename("cavity.off", "cavity.co2.c2v");
    }
};

BOOST_FIXTURE_TEST_SUITE(GePolCavityCO2C2v, GePolCavityCO2C2vTest)

/*! \class GePolCavity
 *  \test \b GePolCavityCO2C2vTest_size tests GePol cavity size for CO2 in C2v symmetry
 */
BOOST_FIXTURE_TEST_CASE(size, GePolCavityCO2C2vTest)
{
    int size = 448;
    int actualSize = cavity.size();
    BOOST_REQUIRE_EQUAL(size, actualSize);
}

/*! \class GePolCavity
 *  \test \b GePolCavityCO2C2vTest_irreducible_size tests GePol cavity irreducible size for CO2 in C2v symmetry
 */
BOOST_FIXTURE_TEST_CASE(irreducible_size, GePolCavityCO2C2vTest)
{
    int size = 112;
    int actualSize = cavity.irreducible_size();
    BOOST_REQUIRE_EQUAL(size, actualSize);
}

/*! \class GePolCavity
 *  \test \b GePolCavityCO2C2vTest_area tests GePol cavity surface area for CO2 in C2v symmetry
 */
BOOST_FIXTURE_TEST_CASE(area, GePolCavityCO2C2vTest)
{
    double area = 250.68176442433020;
    double actualArea = cavity.elementArea().sum();
    BOOST_REQUIRE_CLOSE(area, actualArea, 1.0e-10);
}

/*! \class GePolCavity
 *  \test \b GePolCavityCO2C2vTest_volume tests GePol cavity volume for CO2 in C2v symmetry
 */
BOOST_FIXTURE_TEST_CASE(volume, GePolCavityCO2C2vTest)
{
    double volume = 352.55869984340751;
    Eigen::Matrix3Xd elementCenter = cavity.elementCenter();
    Eigen::Matrix3Xd elementNormal = cavity.elementNormal();
    double actualVolume = 0;
    for ( int i = 0; i < cavity.size(); ++i ) {
        actualVolume += cavity.elementArea(i) * elementCenter.col(i).dot(elementNormal.col(
                            i));
    }
    actualVolume /= 3;
    BOOST_REQUIRE_CLOSE(volume, actualVolume, 1.0e-10);
}

BOOST_AUTO_TEST_SUITE_END()

struct GePolCavityCO2D2hTest {
    GePolCavity cavity;
    GePolCavityCO2D2hTest() { SetUp(); }
    void SetUp() {
        Eigen::Vector3d	C1( 0.0000000000, 0.0000000000, 0.0000000000);
        Eigen::Vector3d	O1( 2.1316110791, 0.0000000000, 0.0000000000);
        Eigen::Vector3d	O2(-2.1316110791, 0.0000000000, 0.0000000000);
        std::vector<Sphere> spheres;
        double radiusC = (1.70 * 1.20) / convertBohrToAngstrom;
        double radiusO = (1.52 * 1.20) / convertBohrToAngstrom;
        Sphere sph1(C1, radiusC);
        Sphere sph2(O1, radiusO);
        Sphere sph3(O2, radiusO);
        spheres.push_back(sph1);
        spheres.push_back(sph2);
        spheres.push_back(sph3);
        double area = 0.2 / convertBohr2ToAngstrom2;
        double probeRadius = 1.385 / convertBohrToAngstrom;
        double minRadius = 0.2 / convertBohrToAngstrom;
        // D2h as generated by Oxy, Oxz and Oyz
        Symmetry pGroup = buildGroup(3, 4, 2, 1);
        cavity = GePolCavity(spheres, area, probeRadius, minRadius, pGroup);
        fs::rename("PEDRA.OUT", "PEDRA.co2.d2h");
        fs::rename("cavity.off", "cavity.co2.d2h");
    }
};

BOOST_FIXTURE_TEST_SUITE(GePolCavityCO2D2h, GePolCavityCO2D2hTest)

/*! \class GePolCavity
 *  \test \b GePolCavityCO2D2hTest_size tests GePol cavity size for CO2 in D2h symmetry
 */
BOOST_FIXTURE_TEST_CASE(size, GePolCavityCO2D2hTest)
{
    int size = 448;
    int actualSize = cavity.size();
    BOOST_REQUIRE_EQUAL(size, actualSize);
}

/*! \class GePolCavity
 *  \test \b GePolCavityCO2D2hTest_irreducible_size tests GePol cavity irreducible size for CO2 in D2h symmetry
 */
BOOST_FIXTURE_TEST_CASE(irreducible_size, GePolCavityCO2D2hTest)
{
    int size = 56;
    int actualSize = cavity.irreducible_size();
    BOOST_REQUIRE_EQUAL(size, actualSize);
}

/*! \class GePolCavity
 *  \test \b GePolCavityCO2D2hTest_area tests GePol cavity surface area for CO2 in D2h symmetry
 */
BOOST_FIXTURE_TEST_CASE(area, GePolCavityCO2D2hTest)
{
    double area = 250.68176442433020;
    double actualArea = cavity.elementArea().sum();
    BOOST_REQUIRE_CLOSE(area, actualArea, 1.0e-10);
}

/*! \class GePolCavity
 *  \test \b GePolCavityCO2D2hTest_volume tests GePol cavity volume for CO2 in D2h symmetry
 */
BOOST_FIXTURE_TEST_CASE(volume, GePolCavityCO2D2hTest)
{
    double volume = 352.55869984340751;
    Eigen::Matrix3Xd elementCenter = cavity.elementCenter();
    Eigen::Matrix3Xd elementNormal = cavity.elementNormal();
    double actualVolume = 0;
    for ( int i = 0; i < cavity.size(); ++i ) {
        actualVolume += cavity.elementArea(i) * elementCenter.col(i).dot(elementNormal.col(
                            i));
    }
    actualVolume /= 3;
    BOOST_REQUIRE_CLOSE(volume, actualVolume, 1.0e-10);
}

BOOST_AUTO_TEST_SUITE_END()
