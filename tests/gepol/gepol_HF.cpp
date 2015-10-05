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

#define BOOST_TEST_MODULE GePolCavityNH3

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include "Config.hpp"

#include <Eigen/Core>

#include "GePolCavity.hpp"
#include "LoggerInterface.hpp"
#include "Molecule.hpp"
#include "PhysicalConstants.hpp"
#include "TestingMolecules.hpp"
#include "Symmetry.hpp"

struct GePolCavityNH3Test {
protected:
    GePolCavity cavity;
    GePolCavityNH3Test() { SetUp(); }
    void SetUp() {
        Molecule molec = HF();
        double area = 0.02 / convertBohr2ToAngstrom2;
        double probeRadius = 0.0 / convertBohrToAngstrom;
        double minRadius = 0.2 / convertBohrToAngstrom;
        cavity = GePolCavity(molec, area, probeRadius, minRadius);
    }
};

BOOST_FIXTURE_TEST_CASE(size, GePolCavityNH3Test)
{
    int ref_size = 1688;
    int size = cavity.size();
    BOOST_REQUIRE_EQUAL(size, ref_size);
}

BOOST_FIXTURE_TEST_CASE(area, GePolCavityNH3Test)
{
    double ref_area = 110.64517236179323;
    double area = cavity.elementArea().sum();
    BOOST_REQUIRE_CLOSE(area, ref_area, 1.0e-09);
}

BOOST_FIXTURE_TEST_CASE(volume, GePolCavityNH3Test)
{
    double ref_volume = 105.98002389543836;
    Eigen::Matrix3Xd elementCenter = cavity.elementCenter();
    Eigen::Matrix3Xd elementNormal = cavity.elementNormal();
    double volume = 0;
    for ( int i = 0; i < cavity.size(); ++i ) {
        volume += cavity.elementArea(i) * elementCenter.col(i).dot(elementNormal.col(i));
    }
    volume /= 3;
    BOOST_REQUIRE_CLOSE(volume, ref_volume, 1.0e-09);
}
