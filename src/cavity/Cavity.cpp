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
#include "MathUtils.hpp"

inline void Cavity::saveCavity(const std::string & fname)
{
	/*
	std::ofstream weights("weights.txt", std::ios_base::out);
	// First line in the weights file is the number of elements.
	// This is for a sanity-check on the save/load operations.
	weights << std::setprecision(std::numeric_limits<double>::digits10) << nElements << std::endl;
	weights << std::setprecision(std::numeric_limits<double>::digits10) << elementArea << std::endl;
	std::ofstream elRadius("element_radius.txt", std::ios_base::out);
	elRadius << std::setprecision(std::numeric_limits<double>::digits10) << elementRadius << std::endl;
	std::ofstream centers("centers.txt", std::ios_base::out);
	centers << std::setprecision(std::numeric_limits<double>::digits10) << elementCenter << std::endl;
	std::ofstream normals("normals.txt", std::ios_base::out);
	normals << std::setprecision(std::numeric_limits<double>::digits10) << elementNormal << std::endl;
	*/
		
	// Write everything in a single .npz binary file
	unsigned int dim = static_cast<unsigned int>(nElements);
	// Write the number of elements, it will be used to check sanity of the save/load operations.
	const unsigned int shape[] = {1};
	cnpy::npz_save(fname, "elements", &dim, shape, 1, "w", false);
	// Write weights
	const unsigned int weights_shape[] = {dim};
	cnpy::npz_save(fname, "weights", elementArea.data(), weights_shape, 1, "a", true);
	// Write element radius
	const unsigned int elRadius_shape[] = {dim};
	cnpy::npz_save(fname, "elRadius", elementRadius.data(), elRadius_shape, 1, "a", true);
	// Write centers
	const unsigned int centers_shape[] = {3, dim};
	cnpy::npz_save(fname, "centers", elementCenter.data(), centers_shape, 2, "a", true);
	// Write normals
	const unsigned int normals_shape[] = {3, dim};
	cnpy::npz_save(fname, "normals", elementNormal.data(), normals_shape, 2, "a", true);
}

inline void Cavity::loadCavity(const std::string & fname)
{
	// If the cavity has been loaded from file, the point group is C1
        pointGroup_ = buildGroup(0);
	// Load the .npz binary file and then traverse it to get the data needed to rebuild the cavity.
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
		elementArea = getFromRawBuffer<double>(dim, 1, raw_weights.data);
	}
	
	// 2. Get the element radius 
        cnpy::NpyArray raw_elRadius = loaded_cavity["elRadius"];
	dim = raw_elRadius.shape[0];
	if (dim != nElements)
	{
		throw std::runtime_error("A problem occurred while loading the cavity. Inconsistent dimension of element radius vector!");
	}
	else
	{
		elementRadius = getFromRawBuffer<double>(dim, 1, raw_elRadius.data);
	}

	// 3. Get the centers
	cnpy::NpyArray raw_centers = loaded_cavity["centers"];
	dim = raw_centers.shape[1];
	if (dim != nElements)
	{
		throw std::runtime_error("A problem occurred while loading the cavity. Inconsistent dimension of centers matrix!");
	}
	else
	{
		elementCenter = getFromRawBuffer<double>(3, dim, raw_centers.data);
	}

	// 4. Get the normal vectors	
	cnpy::NpyArray raw_normals = loaded_cavity["normals"];
	dim = raw_normals.shape[1];
	if (dim != nElements)
	{
		throw std::runtime_error("A problem occurred while loading the cavity. Inconsistent dimension of normals matrix!");
	}
	else
	{
		elementNormal = getFromRawBuffer<double>(3, dim, raw_normals.data);
	}
}
