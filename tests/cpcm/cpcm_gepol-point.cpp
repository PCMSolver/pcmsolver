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
 *     PCMSolver API, see: <http://pcmsolver.readthedocs.org/>
 */
/* pcmsolver_copyright_end */

#include "catch.hpp"

#include <iostream>

#include "Config.hpp"

#include <Eigen/Core>

#include "bi_operators/CollocationIntegrator.hpp"
#include "green/DerivativeTypes.hpp"
#include "cavity/GePolCavity.hpp"
#include "green/Vacuum.hpp"
#include "green/UniformDielectric.hpp"
#include "solver/CPCMSolver.hpp"
#include "utils/TestingMolecules.hpp"

SCENARIO("Test solver for the C-PCM for a point charge and a GePol cavity", "[solver][cpcm][cpcm_gepol-point]")
{
    GIVEN("An isotropic environment modelled as a perfect conductor and a point charge")
    {
        double permittivity = 78.39;
        Vacuum<AD_directional, CollocationIntegrator> gfInside = Vacuum<AD_directional, CollocationIntegrator>();
        UniformDielectric<AD_directional, CollocationIntegrator> gfOutside =
            UniformDielectric<AD_directional, CollocationIntegrator>(permittivity);
        bool symm = true;
        double correction = 0.5;

        double charge = 8.0;
        double totalASC = - charge * (permittivity - 1) / (permittivity + correction);

        /*! \class CPCMSolver
         *  \test \b pointChargeGePol tests CPCMSolver using a point charge with a GePol cavity
         *  The point charge is at the origin.
         */
        WHEN("the point charge is located at the origin")
        {
            Molecule point = dummy<0>(2.929075493);
            double area = 0.4;
            double probeRadius = 0.0;
            double minRadius = 100.0;
            GePolCavity cavity = GePolCavity(point, area, probeRadius, minRadius);
            cavity.saveCavity("point.npz");

            CPCMSolver solver(symm, correction);
            solver.buildSystemMatrix(cavity, gfInside, gfOutside);

            size_t size = cavity.size();
            Eigen::VectorXd fake_mep = computeMEP(cavity.elements(), charge);
            for (size_t i = 0; i < size; ++i) {
                INFO("fake_mep(" << i << ") = " << fake_mep(i));
            }

            THEN("the apparent surface charge is")
            {
                Eigen::VectorXd fake_asc = Eigen::VectorXd::Zero(size);
                fake_asc = solver.computeCharge(fake_mep);
                double totalFakeASC = fake_asc.sum();

                for (size_t i = 0; i < size; ++i) {
                    INFO("fake_asc(" << i << ") = " << fake_asc(i));
                }

                CAPTURE(totalASC);
                CAPTURE(totalFakeASC);
                CAPTURE(totalASC - totalFakeASC);
                REQUIRE(totalASC == Approx(totalFakeASC).epsilon(1.0e-03));
            }
        }

        /*! \class CPCMSolver
         *  \test \b pointChargeShiftedGePol tests CPCMSolver using a point charge with a GePol cavity
         *  The point charge is away from the origin.
         */
        AND_WHEN("the point charge is located away from the origin")
        {
            Eigen::Vector3d origin = 100 * Eigen::Vector3d::Random();
            Molecule point = dummy<0>(2.929075493, origin);
            double area = 0.4;
            double probeRadius = 0.0;
            double minRadius = 100.0;
            GePolCavity cavity = GePolCavity(point, area, probeRadius, minRadius);

            CPCMSolver solver(symm, correction);
            solver.buildSystemMatrix(cavity, gfInside, gfOutside);

            double charge = 8.0;
            size_t size = cavity.size();
            Eigen::VectorXd fake_mep = computeMEP(cavity.elements(), charge, origin);
            for (size_t i = 0; i < size; ++i) {
                INFO("fake_mep(" << i << ") = " << fake_mep(i));
            }

            THEN("the surface charge is")
            {
                Eigen::VectorXd fake_asc = Eigen::VectorXd::Zero(size);
                fake_asc = solver.computeCharge(fake_mep);
                double totalFakeASC = fake_asc.sum();

                for (size_t i = 0; i < size; ++i) {
                    INFO("fake_asc(" << i << ") = " << fake_asc(i));
                }

                CAPTURE(totalASC);
                CAPTURE(totalFakeASC);
                CAPTURE(totalASC - totalFakeASC);
                REQUIRE(totalASC == Approx(totalFakeASC).epsilon(1.0e-03));
            }
        }
    }
}
