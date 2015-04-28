#include <iostream>
#include <ostream>
#include <fstream>

#include <Eigen/Dense>

#include "TanhSphericalDiffuse.hpp"

int main()
{
    int nPoints = 1000;
    double epsInside = 80.0;
    double epsOutside = 2.0;
    Eigen::Vector3d sphereCenter;
    sphereCenter << 0.0, 0.0, 0.0;
    double sphereRadius = 100.0;
    double width = 10.0;

    Eigen::Vector3d source;
    source << 1.0, 0.0, 100.0;
    Eigen::Vector3d probe;
    probe << 0.0, 0.0, 0.0;

    TanhSphericalDiffuse gf(epsInside, epsOutside, width, sphereRadius);

    double zMin = 90.0;
    double zMax = 140.0;
    double step = (zMax - zMin) / nPoints;
    std::ofstream out;
    out.open("gf_plot.log");
    out << "#    Distance             gf_value " << std::endl;
    out.precision(16);
    for (int i = 0; i < nPoints; ++i) {
        probe(2) = zMin + i*step;
        double gf_value = gf.function(source, probe);
        out << probe(2) << "       " << gf_value << std::endl;
    }
    out.close();
}
