#define BOOST_TEST_MODULE CPCMSolverpointChargeGePol

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <iostream>

#include "Config.hpp"

#include "EigenPimpl.hpp"

#include "CollocationIntegrator.hpp"
#include "CPCMSolver.hpp"
#include "DerivativeTypes.hpp"
#include "GePolCavity.hpp"
#include "UniformDielectric.hpp"
#include "Vacuum.hpp"
#include "Symmetry.hpp"

/*! \class CPCMSolver
 *  \test \b pointChargeGePol tests CPCMSolver using a point charge with a GePol cavity
 */
BOOST_AUTO_TEST_CASE(pointChargeGePol)
{
    // Set up cavity
    Eigen::Vector3d N(0.0, 0.0, 0.0);
    std::vector<Sphere> spheres;
    Sphere sph1(N, 2.929075493);
    spheres.push_back(sph1);
    double area = 0.4;
    double probeRadius = 0.0;
    double minRadius = 100.0;
    // C1
    Symmetry group = buildGroup(0, 0, 0, 0);
    GePolCavity cavity = GePolCavity(spheres, area, probeRadius, minRadius, group);
    cavity.saveCavity("point.npz");
    // The point charge is located at the origin.
    // The potential at cavity point s_I is Q/|s_I|
    
    CollocationIntegrator * diag = new CollocationIntegrator();
    double permittivity = 78.39;
    Vacuum<AD_directional> * gfInside = new Vacuum<AD_directional>(diag);
    UniformDielectric<AD_directional> * gfOutside = new
    UniformDielectric<AD_directional>(permittivity, diag);
    bool symm = true;
    double correction = 0.5;
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
    // for CPCM it will be -Q*[(epsilon-1)/epsilon + correction]
    Eigen::VectorXd fake_asc = Eigen::VectorXd::Zero(size);
    solver.compCharge(fake_mep, fake_asc);
    double totalASC = - charge * (permittivity - 1) / (permittivity + correction);
    double totalFakeASC = fake_asc.sum();
    std::cout << "totalASC - totalFakeASC = " << totalASC - totalFakeASC << std::endl;
    BOOST_REQUIRE_CLOSE(totalASC, totalFakeASC, 4e-02);
}
