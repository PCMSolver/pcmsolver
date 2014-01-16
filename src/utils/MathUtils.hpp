#ifndef MATHUTILS_HPP
#define MATHUTILS_HPP

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
}

#endif // MATHUTILS_HPP
