#include "Cavity.hpp"

#include <iomanip>
#include <iostream>
#include <fstream>
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

#include "npy.hpp"

void Cavity::writeCavityFile(bool isNPY)
{
	std::ofstream cavity_file("cavity_spec.txt", std::ios_base::out);
	cavity_file << nElements << std::endl;
	cavity_file << elementArea << std::endl;
	cavity_file << elementCenter << std::endl;
	cavity_file << elementNormal << std::endl;

	if (isNPY)	
	{
		// Write weights
		std::ofstream weights_file ("weights.npy",std::ofstream::binary);
		npy_write_matrix(weights_file, nElements, 1, elementArea.data());
		// Write centers
		std::ofstream centers_file ("centers.npy",std::ofstream::binary);
		npy_write_matrix(centers_file, 3, nElements, elementCenter.data(), true);
		// Write normals
		std::ofstream normals_file ("normals.npy",std::ofstream::binary);
		npy_write_matrix(normals_file, 3, nElements, elementNormal.data(), true);
	}
}

void Cavity::readCavityFile(bool isNPY)
{
	if (!isNPY)
	{
	}
	else
	{
	}
}
