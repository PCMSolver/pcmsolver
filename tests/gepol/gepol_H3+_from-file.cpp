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

#include <vector>
#include <cmath>

#include "Config.hpp"

#include <Eigen/Core>

#include "GePolCavity.hpp"

TEST_CASE("Restart GePol cavity for an H3+ molecule in C1 symmetry", "[gepol][gepol_H3+_from-file]")
{
    GePolCavity cavity;
    cavity.loadCavity("h3+.npz");

    /*! \class GePolCavity
     *  \test \b GePolCavityH3RestartTest_size tests GePol cavity size for H3+ loading the cavity from a .npz file
     */
    SECTION("Test size")
    {
        int size = 312;
        size_t actualSize = cavity.size();
        REQUIRE(size == actualSize);
    }

    /*! \class GePolCavity
     *  \test \b GePolCavityH3RestartTest_area tests GePol cavity surface area for H3+ loading the cavity from from a .npz file
     */
    SECTION("Test surface area")
    {
        double area = 178.74700256125493;
        double actualArea = cavity.elementArea().sum();
        REQUIRE(area == Approx(actualArea));
    }

    /*! \class GePolCavity
     *  \test \b GePolCavityH3RestartTest_volume tests GePol cavity volume for H3+ loading the cavity from from a .npz file
     */
    SECTION("Test volume")
    {
        double volume = 196.4736029455637;
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
