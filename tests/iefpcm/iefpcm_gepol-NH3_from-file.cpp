#include <iostream>

#include "Config.hpp"

#include "EigenPimpl.hpp"

#include "DerivativeTypes.hpp"
#include "GePolCavity.hpp"
#include "Vacuum.hpp"
#include "UniformDielectric.hpp"
#include "IEFSolver.hpp"

#include "gtestPimpl.hpp"

TEST(IEFSolver, NH3GePol)
{
    Eigen::Vector3d N( -0.000000000,   -0.104038047,    0.000000000);
    Eigen::Vector3d H1(-0.901584415,    0.481847022,   -1.561590016);
    Eigen::Vector3d H2(-0.901584415,    0.481847022,    1.561590016);
    Eigen::Vector3d H3( 1.803168833,    0.481847022,    0.000000000);
    // Load  cavity
    GePolCavity cavity;
    cavity.loadCavity("nh3.npz");

    double permittivity = 78.39;
    Vacuum<AD_directional> * gfInside = new Vacuum<AD_directional>();
    UniformDielectric<AD_directional> * gfOutside = new
    UniformDielectric<AD_directional>(permittivity);
    bool symm = true;
    IEFSolver solver(gfInside, gfOutside, symm);
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
    // The total ASC for a dielectric is -Q*[(epsilon-1)/epsilon]
    Eigen::VectorXd fake_asc = Eigen::VectorXd::Zero(size);
    solver.compCharge(fake_mep, fake_asc);
    double totalASC = - (Ncharge + 3.0 * Hcharge) * (permittivity - 1) / permittivity;
    double totalFakeASC = fake_asc.sum();
    std::cout << "totalASC - totalFakeASC = " << totalASC - totalFakeASC << std::endl;
    EXPECT_NEAR(totalASC, totalFakeASC, 3e-3);
}
