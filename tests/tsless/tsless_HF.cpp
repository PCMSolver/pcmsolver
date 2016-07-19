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

#include "catch.hpp"

#include <Eigen/Core>

#include "cavity/TsLessCavity.hpp"
#include "utils/Molecule.hpp"
#include "TestingMolecules.hpp"

TEST_CASE("TsLess cavity for an hydrogen fluoride molecule", "[tsless][tsless_HF]")
{
    Molecule molec = HF();
    double area = 0.02 / bohr2ToAngstrom2();
    double probeRadius = 0.0 / bohrToAngstrom();
    double minRadius = 0.2 / bohrToAngstrom();
    double minDistance = 0.01;
    int derOrder = 4;
    TsLessCavity cavity(molec, area, probeRadius, minRadius, minDistance, derOrder);

    /*! \class TsLessCavity
     *  \test \b TsLessCavityHFTest_size tests TsLess cavity size for ammonia
     */
    SECTION("Test size")
    {
        size_t ref_size = 1498;
        size_t size = cavity.size();
        REQUIRE(size == ref_size);
    }

    /*! \class TsLessCavity
     *  \test \b TsLessCavityHFTest_area tests TsLess cavity surface area for ammonia
     */
    SECTION("Test surface area")
    {
        double ref_area = 111.07794350824591;
        double area = cavity.elementArea().sum();
        REQUIRE(area == Approx(ref_area));
    }

    /*! \class TsLessCavity
     *  \test \b TsLessCavityHFTest_volume tests TsLess cavity volume for ammonia
     */
    SECTION("Test volume")
    {
        double ref_volume = 106.59560252994531;
        Eigen::Matrix3Xd elementCenter = cavity.elementCenter();
        Eigen::Matrix3Xd elementNormal = cavity.elementNormal();
        double volume = 0;
        for ( size_t i = 0; i < cavity.size(); ++i ) {
            volume += cavity.elementArea(i) * elementCenter.col(i).dot(elementNormal.col(i));
        }
        volume /= 3;
        REQUIRE(volume == Approx(ref_volume));
    }
}