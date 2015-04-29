#include <iostream>
#include <ostream>
#include <fstream>

#include <Eigen/Dense>

#include "DerivativeTypes.hpp"
#include "UniformDielectric.hpp"
#include "TanhSphericalDiffuse.hpp"

int main()
{
    int nPoints = 10000;
    double epsInside = 80.0;
    double epsOutside = 1.0;
    Eigen::Vector3d sphereCenter;
    sphereCenter << 0.0, 0.0, 0.0;
    double sphereRadius = 100.0;
    double width = 10.0;

    Eigen::Vector3d source;
    source << 0.5, 0.0, 0.0;
    Eigen::Vector3d probe;
    probe << 50.0, 0.0, 0.0;

    double zMin = 0.0;
    double zMax = 200.0;
    double step = (zMax - zMin) / nPoints;

    TanhSphericalDiffuse gf(epsInside, epsOutside, width, sphereRadius);
    // Plot Green's function, image potential and separation coefficient
    std::ofstream out;
    out.open("gf_spherical.log");
    out << "#" << '\t' << "Distance" << '\t' << "gf_value" << '\t' << "image" << '\t' << "coefficient" << std::endl;
    out.precision(16);
    for (int i = 0; i < nPoints; ++i) {
        probe(2) = zMin + i*step;
        out << '\t' << probe(2) << '\t' << gf.function(source, probe) << '\t' << gf.imagePotential(source, probe) << '\t' << gf.Coulomb(source, probe) << std::endl;
    }
    out.close();

    UniformDielectric<AD_directional> gf_inside(epsInside);
    out.open("gf_uniform_inside.log");
    out << "#" << '\t' << "Distance" << '\t' << "gf_value" << std::endl;
    out.precision(16);
    for (int i = 0; i < nPoints; ++i) {
        probe(2) = zMin + i*step;
        out << '\t' << probe(2) << '\t' << gf_inside.function(source, probe) << std::endl;
    }
    out.close();

    UniformDielectric<AD_directional> gf_outside(epsOutside);
    out.open("gf_uniform_outside.log");
    out << "#" << '\t' << "Distance" << '\t' << "gf_value" << std::endl;
    out.precision(16);
    for (int i = 0; i < nPoints; ++i) {
        probe(2) = zMin + i*step;
        out << '\t' << probe(2) << '\t' << gf_outside.function(source, probe) << std::endl;
    }
    out.close();
}
