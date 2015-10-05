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

#include "catch.hpp"

#include "Config.hpp"

#include <Eigen/Core>

#include "CollocationIntegrator.hpp"
#include "DerivativeTypes.hpp"
#include "Vacuum.hpp"
#include "UniformDielectric.hpp"
#include "IEFSolver.hpp"
#include "TestingMolecules.hpp"
#include "TsLessCavity.hpp"

SCENARIO("Test solver for the IEFPCM for a point charge and a TsLess cavity", "[solver][iefpcm][iefpcm_tsless-point]")
{
    GIVEN("An isotropic environment and a point charge")
    {
        double permittivity = 78.39;
        Vacuum<AD_directional, CollocationIntegrator> gfInside = Vacuum<AD_directional, CollocationIntegrator>();
        UniformDielectric<AD_directional, CollocationIntegrator> gfOutside =
            UniformDielectric<AD_directional, CollocationIntegrator>(permittivity);
        bool symm = true;

        double charge = 8.0;
        double totalASC = - charge * (permittivity - 1) / permittivity;

        /*! \class IEFSolver
         *  \test \b pointChargeTsLess tests IEFSolver using a point charge with a TsLess cavity
         *  The point charge is at the origin.
         */
        WHEN("the point charge is located at the origin")
        {
            Molecule point = dummy<0>(2.929075493);
            double area = 0.4;
            double probeRadius = 0.0;
            double minRadius = 100.0;
            double minDistance = 0.1;
            int derOrder = 4;
            TsLessCavity cavity = TsLessCavity(point, area, probeRadius, minRadius, minDistance, derOrder);
            cavity.saveCavity("tsless_point.npz");

            IEFSolver solver(symm);
            solver.buildSystemMatrix(cavity, gfInside, gfOutside);

            size_t size = cavity.size();
            Eigen::VectorXd fake_mep = computeMEP(cavity.elements(), charge);
            // The total ASC for a dielectric is -Q*[(epsilon-1)/epsilon]
            Eigen::VectorXd fake_asc = Eigen::VectorXd::Zero(size);
            fake_asc = solver.computeCharge(fake_mep);

            for (size_t i = 0; i < size; ++i) {
                INFO("fake_mep(" << i << ") = " << fake_mep(i));
            }
            for (size_t i = 0; i < size; ++i) {
                INFO("fake_asc(" << i << ") = " << fake_asc(i));
            }

            double totalFakeASC = fake_asc.sum();
            THEN("the apparent surface charge is")
            {
                CAPTURE(totalASC);
                CAPTURE(totalFakeASC);
                CAPTURE(totalASC - totalFakeASC);
                REQUIRE(totalASC == Approx(totalFakeASC).epsilon(1.0e-03));
            }
        }

        /*! \class IEFSolver
         *  \test \b pointChargeShiftedTsLess tests IEFSolver using a point charge with a TsLess cavity
         *  The point charge is away from the origin.
         */
        AND_WHEN("the point charge is located away from the origin")
        {
            Eigen::Vector3d origin = 100 * Eigen::Vector3d::Random();
            Molecule point = dummy<0>(2.929075493, origin);
            double area = 0.4;
            double probeRadius = 0.0;
            double minRadius = 100.0;
            double minDistance = 0.1;
            int derOrder = 4;
            TsLessCavity cavity = TsLessCavity(point, area, probeRadius, minRadius, minDistance, derOrder);

            IEFSolver solver(symm);
            solver.buildSystemMatrix(cavity, gfInside, gfOutside);

            size_t size = cavity.size();
            Eigen::VectorXd fake_mep = computeMEP(cavity.elements(), charge, origin);
            // The total ASC for a dielectric is -Q*[(epsilon-1)/epsilon]
            Eigen::VectorXd fake_asc = Eigen::VectorXd::Zero(size);
            fake_asc = solver.computeCharge(fake_mep);

            for (size_t i = 0; i < size; ++i) {
                INFO("fake_mep(" << i << ") = " << fake_mep(i));
            }
            for (size_t i = 0; i < size; ++i) {
                INFO("fake_asc(" << i << ") = " << fake_asc(i));
            }

            double totalFakeASC = fake_asc.sum();
            THEN("the surface charge is")
            {
                CAPTURE(totalASC);
                CAPTURE(totalFakeASC);
                CAPTURE(totalASC - totalFakeASC);
                REQUIRE(totalASC == Approx(totalFakeASC).epsilon(1.0e-03));
            }
        }
    }
}
