#define BOOST_TEST_MODULE CPCMSolverNH3TsLess

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <iostream>

#include "Config.hpp"

#include <Eigen/Dense>

#include "CollocationIntegrator.hpp"
#include "DerivativeTypes.hpp"
#include "TsLessCavity.hpp"
#include "Vacuum.hpp"
#include "UniformDielectric.hpp"
#include "CPCMSolver.hpp"

/*! \class CPCMSolver
 *  \test \b NH3TsLess tests CPCMSolver using ammonia and a TsLess cavity
 */
BOOST_AUTO_TEST_CASE(NH3TsLess)
{
    // Set up cavity
    Eigen::Vector3d N( -0.000000000,   -0.104038047,    0.000000000);
    Eigen::Vector3d H1(-0.901584415,    0.481847022,   -1.561590016);
    Eigen::Vector3d H2(-0.901584415,    0.481847022,    1.561590016);
    Eigen::Vector3d H3( 1.803168833,    0.481847022,    0.000000000);
    std::vector<Sphere> spheres;
    Sphere sph1(N,  2.929075493);
    Sphere sph2(H1, 2.267671349);
    Sphere sph3(H2, 2.267671349);
    Sphere sph4(H3, 2.267671349);
    spheres.push_back(sph1);
    spheres.push_back(sph2);
    spheres.push_back(sph3);
    spheres.push_back(sph4);
    double area = 0.4;
    TsLessCavity cavity(spheres, area);

    CollocationIntegrator * diag = new CollocationIntegrator();
    double permittivity = 78.39;
    Vacuum<AD_directional> * gfInside = new Vacuum<AD_directional>(diag);
    UniformDielectric<AD_directional> * gfOutside = new
    UniformDielectric<AD_directional>(permittivity, diag);
    bool symm = true;
    double correction = 0.0;
    CPCMSolver solver(gfInside, gfOutside, symm, correction);
    solver.buildSystemMatrix(cavity);

    double Ncharge = 7.0;
    double Hcharge = 1.0;
    int size = cavity.size();
    Eigen::VectorXd fake_mep = Eigen::VectorXd::Zero(size);
    for (int i = 0; i < size; ++i) {
        Eigen::Vector3d center = cavity.elementCenter(i);
        double Ndistance = (center - N).norm();
        double H1distance = (center - H1).norm();
        double H2distance = (center - H2).norm();
        double H3distance = (center - H3).norm();
        fake_mep(i) = Ncharge / Ndistance + Hcharge / H1distance + Hcharge / H2distance +
                      Hcharge / H3distance;
    }
    // The total ASC for a conductor is -Q
    // for CPCM it will be -Q*[(epsilon-1)/epsilon}
    Eigen::VectorXd fake_asc = Eigen::VectorXd::Zero(size);
    solver.compCharge(fake_mep, fake_asc);
    double totalASC = - (Ncharge + 3.0 * Hcharge) * (permittivity - 1) / permittivity;
    double totalFakeASC = fake_asc.sum();
    std::cout << "totalASC - totalFakeASC = " << totalASC - totalFakeASC << std::endl;
    BOOST_REQUIRE_CLOSE(totalASC, totalFakeASC, 4e-02);
}
