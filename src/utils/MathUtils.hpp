#ifndef MATHUTILS_HPP
#define MATHUTILS_HPP

#include <iosfwd>
#include <vector>

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

/*! \file MathUtils.hpp
 *  \brief Some math utilities.
 *  \author Roberto Di Remigio
 *  \date 2014
 */

/*! \brief Returns an Eigen matrix of type T, with dimensions _rows*_columns.
 *  \param _rows the number of rows.
 *  \param _columns the number of columns.
 *  \param _inData the raw data buffer.
 *  \tparam T the data type of the matrix to be returned.
 *
 *  Warning! This function assumes that the raw buffer is in column-major order
 *  as in Fortran.
 */

template <typename T>
inline Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> getFromRawBuffer(size_t _rows, size_t _columns, void * _inData)
{
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> _outData;
	_outData.resize(_rows, _columns);
	T * data = reinterpret_cast<T*>(_inData);
	Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> > mapped_data(data, _rows, _columns);
 	_outData = mapped_data;
	return _outData;
};

/*! \brief Packs a block-diagonal matrix.
 *  
 * 
 *  We currently assume that all the blocks are square and have the same dimensionality.
 *  This data type is to be used in conjuction with symmetry handling in the solver.
 */

template <typename T>
class BlockDiagonalMatrix
{
	private:
		typedef std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> > blockVector;
	private:
		blockVector blockedMatrix_;
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> fullMatrix_;
		int nrBlocks_;
		int dimBlock_;
		std::ostream & printObject(std::ostream & os)
		{
                	os << "Matrix has " << nrBlocks_ << " square blocks of dimension " << dimBlock_ << std::endl;
                	for (int i = 0; i < nrBlocks_; ++i)
                	{
                		os << "Block number " << i << std::endl;
                		os << blockedMatrix_[i] << std::endl;
                	}
                	return os;
                }
	public:
		BlockDiagonalMatrix() {}
		/*!
		 *  Packs a block diagonal matrix given the full matrix.
		 */
		BlockDiagonalMatrix(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> & matrix, int nrBlocks, int dimBlock)
			: fullMatrix_(matrix), nrBlocks_(nrBlocks), dimBlock_(dimBlock)
		{	
			// The dimension of the full matrix is (nrBlocks * dimBlock)
			// The full matrix is assumed to be square, with square blocks on the diagonal
			// all with the same dimension.
			int j = 0;
			for (int i = 0; i < nrBlocks_; ++i)
			{
				blockedMatrix_.push_back(fullMatrix_.block(j, j, dimBlock_, dimBlock_));
				j += dimBlock_;
			}
		}
         	friend std::ostream & operator<<(std::ostream & os, BlockDiagonalMatrix<T> & bd)
	 	{
			return bd.printObject(os);
         	}
};

#endif // MATHUTILS_HPP
