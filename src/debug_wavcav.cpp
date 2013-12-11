#include <iostream>
#include <fstream>
#include <string>

#include <Eigen/Dense>

#include "PWCSolver.hpp"
#include "Vacuum.hpp"
#include "UniformDielectric.hpp"
#include "WaveletCavity.hpp"

int main(int argc, char** argv){
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
	double probeRadius = 1.385; // Probe Radius for water
	int patchLevel = 2;
    double coarsity = 0.5;
	WaveletCavity cavity(spheres, probeRadius, patchLevel, coarsity);
	cavity.readCavity("molec_dyadic.dat");
	double permittivity = 78.39;
	Vacuum * gfInside = new Vacuum(2); // Automatic directional derivative
	UniformDielectric * gfOutside = new UniformDielectric(2, permittivity);
	int firstKind = 0;
	PWCSolver solver(gfInside, gfOutside, firstKind);
	solver.buildSystemMatrix(cavity);
	cavity.uploadPoints(solver.getQuadratureLevel(), solver.getT_(), false);
}
