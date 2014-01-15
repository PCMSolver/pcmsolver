#include "Cavity.hpp"

#include <iomanip>
#include <iostream>
#include <fstream>
#include <limits>
#include <stdexcept>

#include "Config.hpp"
#include "FCMangle.hpp"

// Disable obnoxious warnings from Eigen headers
#if defined (__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wall" 
#pragma GCC diagnostic ignored "-Weffc++" 
#pragma GCC diagnostic ignored "-Wextra"
#include <Eigen/Dense>
#pragma GCC diagnostic pop
#elif (__INTEL_COMPILER)
#pragma warning push
#pragma warning disable "-Wall"
#include <Eigen/Dense>
#pragma warning pop
#endif

#include "cnpy.hpp"

void Cavity::saveCavity(const std::string & fname)
{
	std::ofstream weights("weights.txt", std::ios_base::out);
	weights << std::setprecision(std::numeric_limits<double>::digits10) << elementArea << std::endl;
	std::ofstream centers("centers.txt", std::ios_base::out);
	centers << std::setprecision(std::numeric_limits<double>::digits10) << elementCenter << std::endl;
	std::ofstream normals("normals.txt", std::ios_base::out);
	normals << std::setprecision(std::numeric_limits<double>::digits10) << elementNormal << std::endl;
		
	// Write everything in a single .npz binary file
	// Write the number of elements, it will be used to check sanity of the save/load operations.
	const unsigned int shape[] = {1};
	cnpy::npz_save(fname, "elements", &nElements, shape, 1, "w", false);
	// Write weights
	const unsigned int weights_shape[] = {nElements};
	cnpy::npz_save(fname, "weights", elementArea.data(), weights_shape, 1, "a", true);
	// Write centers
	const unsigned int centers_shape[] = {3, nElements};
	cnpy::npz_save(fname, "centers", elementCenter.data(), centers_shape, 2, "a", true);
	// Write normals
	const unsigned int normals_shape[] = {3, nElements};
	cnpy::npz_save(fname, "normals", elementNormal.data(), normals_shape, 2, "a", true);
}

void Cavity::loadCavity(const std::string & fname)
{	
	// Load the .npz binary file and then traverse it to get the data
	// needed to rebuild the cavity.
	cnpy::npz_t loaded_cavity = cnpy::npz_load(fname);
	// 0. Get the number of elements
	cnpy::NpyArray raw_ele = loaded_cavity["elements"];
	int * ne = reinterpret_cast<int*>(raw_ele.data);
	nElements = *ne; 

	// 1. Get the weights
        cnpy::NpyArray raw_weights = loaded_cavity["weights"];
	int dim = raw_weights.shape[0];
	if (dim != nElements)
	{
		throw std::runtime_error("A problem occurred while loading the cavity. Inconsistent dimension of weights vector!");
	}
	else
	{
		elementArea.resize(dim);
		double * loaded_weights = reinterpret_cast<double*>(raw_weights.data);
		Eigen::Map<Eigen::VectorXd> w(loaded_weights, dim, 1);
		elementArea = w;
	}

	// 2. Get the centers
	cnpy::NpyArray raw_centers = loaded_cavity["centers"];
	dim = raw_centers.shape[1];
	if (dim != nElements)
	{
		throw std::runtime_error("A problem occurred while loading the cavity. Inconsistent dimension of centers matrix!");
	}
	else
	{
		elementCenter.resize(Eigen::NoChange, dim);
		double * loaded_centers = reinterpret_cast<double*>(raw_centers.data);
		Eigen::Map<Eigen::Matrix3Xd> c(loaded_centers, 3, dim);
		elementCenter = c;
	}

	// 3. Get the normal vectors	
	cnpy::NpyArray raw_normals = loaded_cavity["normals"];
	dim = raw_normals.shape[1];
	if (dim != nElements)
	{
		throw std::runtime_error("A problem occurred while loading the cavity. Inconsistent dimension of normals matrix!");
	}
	else
	{
		elementNormal.resize(Eigen::NoChange, dim);
		double * loaded_normals = reinterpret_cast<double*>(raw_normals.data);
		Eigen::Map<Eigen::Matrix3Xd> n(loaded_normals, 3, dim);
		elementNormal = n;
	}
}
