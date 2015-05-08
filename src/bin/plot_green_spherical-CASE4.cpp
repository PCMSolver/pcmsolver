#include <iostream>
#include <ostream>
#include <fstream>

#include <Eigen/Dense>

#include "DerivativeTypes.hpp"
#include "UniformDielectric.hpp"
#include "TanhSphericalDiffuse.hpp"
#include "LoggerInterface.hpp"

int main()
{
    int nPoints = 10000;
    double epsInside = 80.0;
    double epsOutside = 2.0;
    Eigen::Vector3d sphereCenter;
    sphereCenter << 0.0, 0.0, 0.0;
    double sphereRadius = 100.0;
    double width = 10.0;

    Eigen::Vector3d source;
    source << 0.0, 0.0, 1.0;
    Eigen::Vector3d probe;
    probe << 0.0, 0.0, 0.0;

    double zMin =  10.0;
    double zMax = 300.0;
    double step = (zMax - zMin) / nPoints;

    TanhSphericalDiffuse gf(epsInside, epsOutside, width, sphereRadius);
    LOG(gf);
    std::ofstream out;
    out.open("gf_spherical.log");
    out << "#" << '\t' << "Distance" << '\t' << "gf_value" << '\t' << "image" << '\t' << "coefficient" << '\t' << "derivative" << std::endl;
    out.precision(16);
    for (int i = 0; i < nPoints; ++i) {
        probe(2) = zMin + i*step;
        Eigen::Vector3d direction = probe;
        direction.normalize();
        out << '\t' << probe(2)
            << '\t' << gf.function(source, probe)
            << '\t' << gf.imagePotential(source, probe)
            << '\t' << gf.Coulomb(source, probe)
            << '\t' << gf.derivative(direction, source, probe)
            << '\t' << gf.CoulombDerivative(direction, source, probe)
            << '\t' << gf.imagePotentialDerivative(direction, source, probe) << std::endl;
    }
    out.close();
    LOG_TIME;

    UniformDielectric<AD_directional> gf_inside(epsInside);
    out.open("gf_uniform_inside.log");
    out << "#" << '\t' << "Distance" << '\t' << "gf_value" << std::endl;
    out.precision(16);
    for (int i = 0; i < nPoints; ++i) {
        probe(2) = zMin + i*step;
        Eigen::Vector3d direction = probe;
        direction.normalize();
        out << '\t' << probe(2)
            << '\t' << gf_inside.function(source, probe)
            << '\t' << gf_inside.derivative(direction, source, probe) << std::endl;
    }
    out.close();

    UniformDielectric<AD_directional> gf_outside(epsOutside);
    out.open("gf_uniform_outside.log");
    out << "#" << '\t' << "Distance" << '\t' << "gf_value" << std::endl;
    out.precision(16);
    for (int i = 0; i < nPoints; ++i) {
        probe(2) = zMin + i*step;
        Eigen::Vector3d direction = probe;
        direction.normalize();
        out << '\t' << probe(2)
            << '\t' << gf_outside.function(source, probe)
            << '\t' << gf_outside.derivative(direction, source, probe) << std::endl;
    }
    out.close();
}
