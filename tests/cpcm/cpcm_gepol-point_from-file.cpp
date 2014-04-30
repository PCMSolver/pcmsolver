#define BOOST_TEST_MODULE CPCMSolverpointChargeGePolRestart

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <iostream>

#include "Config.hpp"

#include "EigenPimpl.hpp"

#include "DerivativeTypes.hpp"
#include "GePolCavity.hpp"
#include "Vacuum.hpp"
#include "UniformDielectric.hpp"
#include "CPCMSolver.hpp"

/*! \class CPCMSolver
 *  \test \b pointChargeGePolRestart tests CPCMSolver using a point charge with a GePol cavity read from .npz file
 */
BOOST_AUTO_TEST_CASE(pointChargeGePolRestart)
{
    // Set up cavity
    GePolCavity cavity;
    cavity.loadCavity("point.npz");
    // The point charge is located at the origin.
    // The potential at cavity point s_I is Q/|s_I|
    double permittivity = 78.39;
    Vacuum<AD_directional> * gfInside = new Vacuum<AD_directional>();
    UniformDielectric<AD_directional> * gfOutside = new
    UniformDielectric<AD_directional>(permittivity);
    bool symm = true;
    double correction = 0.0;
    CPCMSolver solver(gfInside, gfOutside, symm, correction);
    solver.buildSystemMatrix(cavity);

    double charge = 8.0;
    int size = cavity.size();
    Eigen::VectorXd fake_mep = Eigen::VectorXd::Zero(size);
    for (int i = 0; i < size; ++i) {
        Eigen::Vector3d center = cavity.elementCenter(i);
        double distance = center.norm();
        fake_mep(i) = charge / distance;
    }
    // The total ASC for a conductor is -Q
    // for CPCM it will be -Q*[(epsilon-1)/epsilon}
    Eigen::VectorXd fake_asc = Eigen::VectorXd::Zero(size);
    solver.compCharge(fake_mep, fake_asc);
    double totalASC = - charge * (permittivity - 1) / permittivity;
    double totalFakeASC = fake_asc.sum();
    std::cout << "totalASC - totalFakeASC = " << totalASC - totalFakeASC << std::endl;
    BOOST_REQUIRE_CLOSE(totalASC, totalFakeASC, 4e-02);
}
