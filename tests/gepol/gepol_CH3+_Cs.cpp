/* pcmsolver_copyright_start */
/*
 *     PCMSolver, an API for the Polarizable Continuum Model
 *     Copyright (C) 2013-2015 Roberto Di Remigio, Luca Frediani and contributors
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

#include <catch.hpp>

#include <cmath>

#include "Config.hpp"

#include <Eigen/Core>

#include "GePolCavity.hpp"
#include "Molecule.hpp"
#include "PhysicalConstants.hpp"
#include "Symmetry.hpp"
#include "TestingMolecules.hpp"

TEST_CASE("GePol cavity for the CH3+ molecule in Cs symmetry", "[gepol][gepol_CH3+_Cs]")
{
    double area = 0.2 / convertBohr2ToAngstrom2;
    double probeRadius = 1.385 / convertBohrToAngstrom;
    // Addition of spheres is enabled, but will not happen in this particular case
    double minRadius = 0.2 / convertBohrToAngstrom;
    Molecule molec = CH3();
    GePolCavity cavity = GePolCavity(molec, area, probeRadius, minRadius, "cs");

    /*! \class GePolCavity
     *  \test \b GePolCavityCsAddTest_size tests GePol cavity size for CH3+ in Cs symmetry with added spheres
     */
    SECTION("Test size")
    {
        int size = 384;
        size_t actualSize = cavity.size();
        REQUIRE(size == actualSize);
    }

    /*! \class GePolCavity
     *  \test \b GePolCavityCsAddTest_irreducible_size tests GePol cavity irreducible size for CH3+ in Cs symmetry with added spheres
     */
    SECTION("Test irreducible size")
    {
        int size = 192;
        int actualSize = cavity.irreducible_size();
        REQUIRE(size == actualSize);
    }

    /*! \class GePolCavity
     *  \test \b GePolCavityCsAddTest_area tests GePol cavity surface area for CH3+ in Cs symmetry with added spheres
     */
    SECTION("Test surface area")
    {
        double area = 211.86178059383573;
        double actualArea = cavity.elementArea().sum();
        REQUIRE(area == Approx(actualArea));
    }

    /*! \class GePolCavity
     *  \test \b GePolCavityCsAddTest_volume tests GePol cavity volume for CH3+ in Cs symmetry with added spheres
     */
    SECTION("Test volume")
    {
        double volume = 278.95706420724309;
        Eigen::Matrix3Xd elementCenter = cavity.elementCenter();
        Eigen::Matrix3Xd elementNormal = cavity.elementNormal();
        double actualVolume = 0;
        for ( size_t i = 0; i < cavity.size(); ++i ) {
            actualVolume += cavity.elementArea(i) * elementCenter.col(i).dot(elementNormal.col(
                        i));
        }
        actualVolume /= 3;
        REQUIRE(volume == Approx(actualVolume));
    }
}
