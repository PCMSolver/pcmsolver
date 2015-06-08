#include <iostream>
#include <ostream>
#include <fstream>
#include <sstream>
#include <string>

#include <Eigen/Dense>

#include "CollocationIntegrator.hpp"
#include "GeneralPurpose.hpp"
#include "LoggerInterface.hpp"
#include "OneLayerTanh.hpp"
#include "SphericalDiffuse.hpp"
#include "UniformDielectric.hpp"

std::vector<double> generate_grid(double, double, size_t);
std::string header();
void case1();
void case2();
void case3();
void case4();

int main()
{
    case1();
//    case2();
//    case3();
//    case4();

    return EXIT_SUCCESS;
}

std::vector<double> generate_grid(double min_val, double max_val, size_t nPoints)
{
    double step = (max_val - min_val) / nPoints;
    std::vector<double> grid;
    for (size_t i = 0; i < nPoints; ++i) {
        grid.push_back(min_val + i * step);
    }
    return grid;
}

std::string header()
{
    std::stringstream tmp;
    tmp << "#" << '\t' << "x_2" << '\t' << "z_2" << '\t';
    tmp << "gf_value" << '\t' << "image" << '\t' << "Coulomb" << '\t';
    tmp << "derivativeProbe" << '\t' << "der_image" << '\t' << "der_Coulomb" << std::endl;
    return tmp.str();
}

void case1()
{
    size_t nPoints = 100;
    double epsInside = 80.0;
    double epsOutside = 2.0;
    Eigen::Vector3d sphereCenter = Eigen::Vector3d::Zero();
    double sphereRadius = 100.0;
    double width = 10.0;

    Eigen::Vector3d source;
    source << 4.0, 0.0, 85.0;
    Eigen::Vector3d probe = Eigen::Vector3d::Zero();

    double x_min = 10.0;
    double x_max = 300.0;
    std::vector<double> x = generate_grid(x_min, x_max, nPoints);
    double z_min = 10.0;
    double z_max = 300.0;
    std::vector<double> z = generate_grid(z_min, z_max, nPoints);

    std::ofstream out;

    SphericalDiffuse<CollocationIntegrator, OneLayerTanh> gf(epsInside, epsOutside, width, sphereRadius, sphereCenter);
    LOG(gf);

    out.open("gf_spherical_CASE1.dat");
    out << header();
    out.precision(16);
    for (size_t i = 0; i < nPoints; ++i) {
        for (size_t j = 0; j < nPoints; ++j) {
            probe << x[i], 0.0, z[j];
            out << '\t' << x[i]
                << '\t' << z[j]
                << '\t' << gf.kernelS(source, probe)
                << '\t' << gf.imagePotential(source, probe)
                << '\t' << gf.Coulomb(source, probe)
                << '\t' << gf.derivativeProbe(Eigen::Vector3d::UnitZ(), source, probe)
                << '\t' << gf.CoulombDerivative(Eigen::Vector3d::UnitZ(), source, probe)
                << '\t' << gf.imagePotentialDerivative(Eigen::Vector3d::UnitZ(), source, probe) << std::endl;
        }
    }
    out.close();
    LOG_TIME;

    UniformDielectric<AD_directional, CollocationIntegrator> gf_inside(epsInside);
    out.open("gf_uniform_inside_CASE1.dat");
    out << "#" << '\t' << "x_2" << '\t' << "z_2" << '\t';
    out << "gf_value" << '\t' << "derivativeProbe" << std::endl;
    for (size_t i = 0; i < nPoints; ++i) {
        for (size_t j = 0; j < nPoints; ++j) {
            probe << x[i], 0.0, z[j];
            out << '\t' << x[i]
                << '\t' << z[j]
                << '\t' << gf_inside.kernelS(source, probe)
                << '\t' << gf_inside.derivativeProbe(Eigen::Vector3d::UnitZ(), source, probe) << std::endl;
        }
    }
    out.precision(16);
    out.close();

    UniformDielectric<AD_directional, CollocationIntegrator> gf_outside(epsOutside);
    out.open("gf_uniform_outside_CASE1.dat");
    out << "#" << '\t' << "x_2" << '\t' << "z_2" << '\t';
    out << "gf_value" << '\t' << "derivativeProbe" << std::endl;
    out.precision(16);
    for (size_t i = 0; i < nPoints; ++i) {
        for (size_t j = 0; j < nPoints; ++j) {
            probe << x[i], 0.0, z[j];
            out << '\t' << x[i]
                << '\t' << z[j]
                << '\t' << gf_outside.kernelS(source, probe)
                << '\t' << gf_outside.derivativeProbe(Eigen::Vector3d::UnitZ(), source, probe) << std::endl;
        }
    }
    out.close();
}
