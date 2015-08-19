#include <iostream>
#include <ostream>
#include <fstream>
#include <sstream>
#include <string>

#include <Eigen/Core>

#include "CollocationIntegrator.hpp"
#include "DerivativeTypes.hpp"
#include "GeneralPurpose.hpp"
#include "LoggerInterface.hpp"
#include "OneLayerTanh.hpp"
#include "PhysicalConstants.hpp"
#include "SphericalDiffuse.hpp"
#include "UniformDielectric.hpp"

/*! \brief Generate uniform grid of points in interval
 *  \param[in] min_val lower bound of interval
 *  \param[in] max_val upper bound of interval
 *  \param[in] nPoints number of points in interval
 */
std::vector<double> generate_grid(double, double, size_t);

/*! \brief Plot Green's function
 *  The source charge is at (0.5, 0.5, 0.0) i.e. inside the dielectric sphere
 */
void case1(double, double, double, double, const Eigen::Vector3d &, int, const std::vector<double> &);

/*! \brief Plot Green's function
 *  The source charge is at (50.0, 50.0, 0.0) i.e. outside the dielectric sphere
 */
void case2(double, double, double, double, const Eigen::Vector3d &, int, const std::vector<double> &);

/*! \brief Plot Green's function
 *  The source charge is at (12.0, 12.0, 0.0) i.e. in the transition layer
 */
void case3(double, double, double, double, const Eigen::Vector3d &, int, const std::vector<double> &);

int main()
{
    size_t nPoints = 200;
    double epsInside = 78.39;
    double epsOutside = 1.0;
    double angToBohr = angstromToBohr(2010);
    Eigen::Vector3d sphereCenter = Eigen::Vector3d::Zero();
    sphereCenter << 0.0, 0.0, 0.0;
    sphereCenter *= angToBohr;
    double sphereRadius = 10.0 * angToBohr;
    double width = 5.0 * angToBohr;
    int maxL = 150;

    // Conversion is done later, so that x, y values are printed in Angstrom
    double min = 0.0;
    double max = 100.0;
    std::vector<double> grid = generate_grid(min, max, nPoints);

    case1(epsInside, epsOutside, width, sphereRadius, sphereCenter, maxL, grid);
    case2(epsInside, epsOutside, width, sphereRadius, sphereCenter, maxL, grid);
    case3(epsInside, epsOutside, width, sphereRadius, sphereCenter, maxL, grid);

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

void case1(double epsInside, double epsOutside, double width, double sphereRadius, const Eigen::Vector3d & sphereCenter, int maxL, const std::vector<double> & grid)
{
    double angToBohr = angstromToBohr(2010);
    Eigen::Vector3d source = Eigen::Vector3d::Zero();
    source << 0.5, 0.5, 0.0;
    source *= angToBohr;
    Eigen::Vector3d probe = Eigen::Vector3d::Zero();

    double delta = 0.1;

    std::ofstream out;

    SphericalDiffuse<CollocationIntegrator, OneLayerTanh> gf(epsInside, epsOutside, width, sphereRadius, sphereCenter, maxL);
    LOG(gf);
    gf.toFile("CASE1");
    UniformDielectric<AD_directional, CollocationIntegrator> gf_i(epsInside);
    UniformDielectric<AD_directional, CollocationIntegrator> gf_o(epsOutside);

    out.open("gf_spherical_CASE1.dat");
    out.precision(16);
    for (size_t i = 0; i < grid.size(); ++i) {
        for (size_t j = 0; j < grid.size(); ++j) {
            probe << (grid[i] + delta), (grid[j] + delta), delta;
	    probe *= angToBohr;
            out << '\t' << i
		<< '\t' << j
		<< '\t' << grid[i] + delta
                << '\t' << grid[j] + delta
                << '\t' << gf.Coulomb(source, probe)
                << '\t' << gf.imagePotential(source, probe)
                << '\t' << gf.kernelS(source, probe)
		<< '\t' << (1.0 / gf.coefficientCoulomb(source, probe)) 
		<< '\t' << gf_i.kernelS(source, probe) 
		<< '\t' << gf_o.kernelS(source, probe) 
		<< std::endl;
        }
	out << std::endl;
    }
    out.close();
    LOG_TIME;
}

void case2(double epsInside, double epsOutside, double width, double sphereRadius, const Eigen::Vector3d & sphereCenter, int maxL, const std::vector<double> & grid)
{
    double angToBohr = angstromToBohr(2010);

    Eigen::Vector3d source = Eigen::Vector3d::Zero();
    source << 50.0, 50.0, 0.0;
    source *= angToBohr;
    Eigen::Vector3d probe = Eigen::Vector3d::Zero();

    double delta = 0.1;

    std::ofstream out;

    SphericalDiffuse<CollocationIntegrator, OneLayerTanh> gf(epsInside, epsOutside, width, sphereRadius, sphereCenter, maxL);
    LOG(gf);
    gf.toFile("CASE2");
    UniformDielectric<AD_directional, CollocationIntegrator> gf_i(epsInside);
    UniformDielectric<AD_directional, CollocationIntegrator> gf_o(epsOutside);

    out.open("gf_spherical_CASE2.dat");
    out.precision(16);
    for (size_t i = 0; i < grid.size(); ++i) {
        for (size_t j = 0; j < grid.size(); ++j) {
            probe << (grid[i] + delta), (grid[j] + delta), delta;
	    probe *= angToBohr;
            out << '\t' << i
		<< '\t' << j
		<< '\t' << grid[i] + delta
                << '\t' << grid[j] + delta
                << '\t' << gf.Coulomb(source, probe)
                << '\t' << gf.imagePotential(source, probe)
                << '\t' << gf.kernelS(source, probe)
		<< '\t' << (1.0 / gf.coefficientCoulomb(source, probe)) 
		<< '\t' << gf_i.kernelS(source, probe) 
		<< '\t' << gf_o.kernelS(source, probe) 
		<< std::endl;
        }
	out << std::endl;
    }
    out.close();
    LOG_TIME;
}

void case3(double epsInside, double epsOutside, double width, double sphereRadius, const Eigen::Vector3d & sphereCenter, int maxL, const std::vector<double> & grid)
{
    double angToBohr = angstromToBohr(2010);

    Eigen::Vector3d source = Eigen::Vector3d::Zero();
    source << 12.0, 12.0, 0.0;
    source *= angToBohr;
    Eigen::Vector3d probe = Eigen::Vector3d::Zero();

    double delta = 0.1;

    std::ofstream out;

    SphericalDiffuse<CollocationIntegrator, OneLayerTanh> gf(epsInside, epsOutside, width, sphereRadius, sphereCenter, maxL);
    LOG(gf);
    gf.toFile("CASE3");
    UniformDielectric<AD_directional, CollocationIntegrator> gf_i(epsInside);
    UniformDielectric<AD_directional, CollocationIntegrator> gf_o(epsOutside);

    out.open("gf_spherical_CASE3.dat");
    out.precision(16);
    for (size_t i = 0; i < grid.size(); ++i) {
        for (size_t j = 0; j < grid.size(); ++j) {
            probe << (grid[i] + delta), (grid[j] + delta), delta;
	    probe *= angToBohr;
            out << '\t' << i
		<< '\t' << j
		<< '\t' << grid[i] + delta
                << '\t' << grid[j] + delta
                << '\t' << gf.Coulomb(source, probe)
                << '\t' << gf.imagePotential(source, probe)
                << '\t' << gf.kernelS(source, probe)
		<< '\t' << (1.0 / gf.coefficientCoulomb(source, probe)) 
		<< '\t' << gf_i.kernelS(source, probe) 
		<< '\t' << gf_o.kernelS(source, probe) 
		<< std::endl;
        }
	out << std::endl;
    }
    out.close();
    LOG_TIME;
}
