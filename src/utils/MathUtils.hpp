#ifndef MATHUTILS_HPP
#define MATHUTILS_HPP

#include <cmath>
#include <iostream>

#include "Config.hpp"

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

#include "Symmetry.hpp"

/*! \fn inline void symmetryBlocking(Eigen::MatrixXd & matrix, int cavitySize, int ntsirr, int nr_irrep)
 *  \param[out] matrix the matrix to be block-diagonalized
 *  \param[in]  cavitySize the size of the cavity (size of the matrix)
 *  \param[in]  ntsirr     the size of the irreducible portion of the cavity (size of the blocks)
 *  \param[in]  nr_irrep   the number of irreducible representations (number of blocks)
 */
inline void symmetryBlocking(Eigen::MatrixXd & matrix, int cavitySize, int ntsirr, int nr_irrep)
{
	// This function implements the simmetry-blocking of the PCM
	// matrix due to point group symmetry as reported in:
	// L. Frediani, R. Cammi, C. S. Pomelli, J. Tomasi and K. Ruud, J. Comput.Chem. 25, 375 (2003)
	// u is the character table for the group (t in the paper)
	Eigen::MatrixXd u = Eigen::MatrixXd::Zero(nr_irrep, nr_irrep);
	for (int i = 0; i < nr_irrep; ++i)
	{
		for (int j = 0; j < nr_irrep; ++j)
		{
			u(i, j) = Symmetry::parity(i&j);
		}
	}
	// Naming of indices:
	//     a, b, c, d   run over the total size of the cavity (N)
	//     i, j, k, l   run over the number of irreps (n)
	//     p, q, r, s   run over the irreducible size of the cavity (N/n)
	// Instead of forming U (T in the paper) and then perform the dense
	// multiplication, we multiply block-by-block using just the u matrix.
	//      matrix = U * matrix * Ut; U * Ut = Ut * U = id
	// First half-transformation, i.e. first_half = matrix * Ut
	Eigen::MatrixXd first_half = Eigen::MatrixXd::Zero(cavitySize, cavitySize);
	for (int i = 0; i < nr_irrep; ++i)
	{
		int ioff = i * ntsirr;
		for (int k = 0; k < nr_irrep; ++k)
		{
			int koff = k * ntsirr;
			for (int j = 0; j < nr_irrep; ++j)
			{
				int joff = j * ntsirr;
				double ujk = u(j, k) / nr_irrep;
				for (int p = 0; p < ntsirr; ++p)
				{
					int a = ioff + p;
					for (int q = 0; q < ntsirr; ++q)
					{
						int b = joff + q;
						int c = koff + q;
						first_half(a, c) += matrix(a, b) * ujk;
					}
				}
			}
		}
	}
	// Second half-transformation, i.e. matrix = U * first_half
	matrix.setZero(cavitySize, cavitySize);
	for (int i = 0; i < nr_irrep; ++i)
	{
		int ioff = i * ntsirr;
		for (int k = 0; k < nr_irrep; ++k) 
		{
			int koff = k * ntsirr;
			for (int j = 0; j < nr_irrep; ++j)
			{
				int joff = j * ntsirr;
				double uij = u(i, j);
				for (int p = 0; p < ntsirr; ++p)
				{
					int a = ioff + p;
					int b = joff + p;
					for (int q = 0; q < ntsirr; ++q)
					{
						int c = koff + q;
						matrix(a, c) += uij * first_half(b, c);
					}
				}
			}
		}
	}
	// Traverse the matrix and discard numerical zeros
	for (int a = 0; a < cavitySize; ++a)
	{
		for (int b = 0; b < cavitySize; ++b)
		{
			if ( std::abs(matrix(a, b)) < 1.0e-14 )
			{
				matrix(a, b) = 0.0;	        				
			}
		}
	}
}

/*! \fn inline void symmetryPacking(std::vector<Eigen::MatrixXd> & blockedMatrix, const Eigen::MatrixXd & fullMatrix, int nrBlocks, int dimBlock)
 *  \param[out] blockedMatrix the result of packing fullMatrix
 *  \param[in]  fullMatrix the matrix to be packed
 *  \param[in]  dimBlock the dimension of the square blocks
 *  \param[in]  nrBlocks the number of square blocks
 */
inline void symmetryPacking(std::vector<Eigen::MatrixXd> & blockedMatrix, const Eigen::MatrixXd & fullMatrix, int dimBlock, int nrBlocks)
{
	// This function packs the square block diagonal fullMatrix with nrBlocks of dimension dimBlock 
	// into a std::vector<Eigen::MatrixXd> containing nrBlocks square matrices of dimenion dimBlock.
	int j = 0;
        for (int i = 0; i < nrBlocks; ++i)
	{
		blockedMatrix.push_back(fullMatrix.block(j, j, dimBlock, dimBlock));
		j += dimBlock;
	}
}

#endif // MATHUTILS_HPP
