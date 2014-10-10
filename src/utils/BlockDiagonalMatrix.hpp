/* pcmsolver_copyright_start */
/*
 *     PCMSolver, an API for the Polarizable Continuum Model
 *     Copyright (C) 2013 Roberto Di Remigio, Luca Frediani and contributors
 *
 *     This file is part of PCMSolver.
 *
 *     PCMSolver is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU Lesser General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     PCMSolver is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU Lesser General Public License for more details.
 *
 *     You should have received a copy of the GNU Lesser General Public License
 *     along with PCMSolver.  If not, see <http://www.gnu.org/licenses/>.
 *
 *     For information on the complete list of contributors to the
 *     PCMSolver API, see: <http://pcmsolver.github.io/pcmsolver-doc>
 */
/* pcmsolver_copyright_end */

#ifndef BLOCKDIAGONALMATRIX_HPP
#define BLOCKDIAGONALMATRIX_HPP

#include <iosfwd>
#include <stdexcept>
#include <vector>

#include "Config.hpp"

#include <Eigen/Dense>

/*! \file BlockDiagonalMatrix.hpp
 *  \class BlockDiagonalMatrix
 *  \brief Packs a block-diagonal matrix.
 *  \author Roberto Di Remigio
 *  \date 2014
 *  \tparam         T the data type
 *  \tparam nrBlocks_ the number of square blocks in the matrix
 *  \tparam dimBlock_ the size of the square blocks
 *
 *  We currently assume that all the blocks are square and have the same dimensionality.
 *  There is no check on these assumptions in any of the public methods: any other
 *  use of this class will lead to runtime errors or "slicing" of the original matrix.
 *  This data type is to be used in conjuction with symmetry handling in the solver.
 */

template <typename T, int nrBlocks_, int dimBlock_>
class BlockDiagonalMatrix
{
private:
    typedef Eigen::Matrix<T, dimBlock_, dimBlock_, Eigen::ColMajor> EigenBlock;
    typedef Eigen::Matrix<T, (nrBlocks_*dimBlock_), (nrBlocks_*dimBlock_), Eigen::ColMajor>
    EigenFull;
    typedef std::vector<EigenBlock> blockVector;
private:
    blockVector blockedMatrix_;
    EigenFull fullMatrix_;
    int fullDim_;
    std::ostream & printObject(std::ostream & os) {
        os << "Matrix has " << nrBlocks_ << " square blocks of dimension " << dimBlock_ <<
           std::endl;
        for (int i = 0; i < nrBlocks_; ++i) {
            os << "Block number " << i << std::endl;
            os << blockedMatrix_[i] << std::endl;
        }
        return os;
    }
public:
    BlockDiagonalMatrix() {}
    /*!
     *  Packs a square block diagonal matrix given the full matrix
     *  The dimension of the full matrix is (nrBlocks * dimBlock)
     *  and it is assumed to be square, with square blocks on the diagonal
     *  all with the same dimension.
     */
    BlockDiagonalMatrix(const EigenFull & matrix) : fullMatrix_(matrix),
        fullDim_(nrBlocks_*dimBlock_) {
        int j = 0;
        for (int i = 0; i < nrBlocks_; ++i) {
            blockedMatrix_.push_back(fullMatrix_.block(j, j, dimBlock_, dimBlock_));
            j += dimBlock_;
        }
    }
    /// Copy constructor
    BlockDiagonalMatrix(const BlockDiagonalMatrix & other) : blockedMatrix_
        (other.blockedMatrix_), fullMatrix_(other.fullMatrix_), fullDim_(other.fullDim_) {}
    friend inline void swap(BlockDiagonalMatrix & left, BlockDiagonalMatrix & right);
    inline void swap(BlockDiagonalMatrix & other);
    /// Addition-assignment operator
    BlockDiagonalMatrix & operator+=(const BlockDiagonalMatrix & other);
    /// Subtraction-assignment operator
    BlockDiagonalMatrix & operator-=(const BlockDiagonalMatrix & other);
    /// Multiplication-assignment operator: BlockDiagonalMatrix * BlockDiagonalMatrix case
    BlockDiagonalMatrix & operator*=(const BlockDiagonalMatrix & other);
    /// Multiplication-assignment operator: uniform scaling case
    BlockDiagonalMatrix & operator*=(double scaling);

    ~BlockDiagonalMatrix() {}
    friend std::ostream & operator<<(std::ostream & os,
                                     BlockDiagonalMatrix<T, nrBlocks_, dimBlock_> & bd) {
        return bd.printObject(os);
    }

    int nrBlocks() const { return nrBlocks_ ; }
    int dimBlock() const { return dimBlock_; }
    int fullDim() const { return fullDim_; }
    const EigenFull & fullMatrix() const { return fullMatrix_; }
    void fullMatrix(const EigenFull & matrix) const { matrix = fullMatrix_; }
    const EigenBlock & block(int i) const { return blockedMatrix_[i]; }
    void block(const EigenBlock & matrix, int i) const { matrix = blockedMatrix_[i]; }
};

template <typename T, int nrBlocks_, int dimBlock_>
BlockDiagonalMatrix<T, nrBlocks_, dimBlock_> &
BlockDiagonalMatrix<T, nrBlocks_, dimBlock_>::operator+=(const
        BlockDiagonalMatrix<T, nrBlocks_, dimBlock_> & other)
{
    for (int i = 0; i < nrBlocks_; ++i) {
        blockedMatrix_[i] += other.blockedMatrix_[i];
    }
    return *this;
}

template <typename T, int nrBlocks_, int dimBlock_>
BlockDiagonalMatrix<T, nrBlocks_, dimBlock_> &
BlockDiagonalMatrix<T, nrBlocks_, dimBlock_>::operator-=(const
        BlockDiagonalMatrix<T, nrBlocks_, dimBlock_> & other)
{
    for (int i = 0; i < nrBlocks_; ++i) {
        blockedMatrix_[i] -= other.blockedMatrix_[i];
    }
    return *this;
}

template <typename T, int nrBlocks_, int dimBlock_>
BlockDiagonalMatrix<T, nrBlocks_, dimBlock_> &
BlockDiagonalMatrix<T, nrBlocks_, dimBlock_>::operator*=(const
        BlockDiagonalMatrix<T, nrBlocks_, dimBlock_> & other)
{
    for (int i = 0; i < nrBlocks_; ++i) {
        blockedMatrix_[i] *= other.blockedMatrix_[i];
    }
    return *this;
}

template <typename T, int nrBlocks_, int dimBlock_>
BlockDiagonalMatrix<T, nrBlocks_, dimBlock_> &
BlockDiagonalMatrix<T, nrBlocks_, dimBlock_>::operator*=(double scaling)
{
    for (int i = 0; i < nrBlocks_; ++i) {
        blockedMatrix_[i] *= scaling;
    }
    return *this;
}

/*!
 * \fn inline BlockDiagonalMatrix operator+(BlockDiagonalMatrix left, const BlockDiagonalMatrix & right)
 * \brief Addition operator
 * \param left the left hand side of the addition
 * \param right the right hand side of the addition
 */
template <typename T, int nrBlocks_, int dimBlock_>
inline BlockDiagonalMatrix<T, nrBlocks_, dimBlock_> operator+
(BlockDiagonalMatrix<T, nrBlocks_, dimBlock_> left,
 const BlockDiagonalMatrix<T, nrBlocks_, dimBlock_> & right)
{
    left += right;
    return left;
}

/*!
 * \fn inline BlockDiagonalMatrix operator-(BlockDiagonalMatrix left, const BlockDiagonalMatrix & right)
 * \brief Subtraction operator
 * \param left the left hand side of the subtraction
 * \param right the right hand side of the subtraction
 */
template <typename T, int nrBlocks_, int dimBlock_>
inline BlockDiagonalMatrix<T, nrBlocks_, dimBlock_> operator-
(BlockDiagonalMatrix<T, nrBlocks_, dimBlock_> left,
 const BlockDiagonalMatrix<T, nrBlocks_, dimBlock_> & right)
{
    left -= right;
    return left;
}

/*!
 * \fn inline BlockDiagonalMatrix operator*(BlockDiagonalMatrix left, const BlockDiagonalMatrix & right)
 * \brief Multiplication operator: BlockDiagonalMatrix * BlockDiagonalMatrix version
 * \param left the left hand side of the multiplication
 * \param right the right hand side of the multiplication
 */
template <typename T, int nrBlocks_, int dimBlock_>
inline BlockDiagonalMatrix<T, nrBlocks_, dimBlock_> operator*
(BlockDiagonalMatrix<T, nrBlocks_, dimBlock_> left,
 const BlockDiagonalMatrix<T, nrBlocks_, dimBlock_> & right)
{
    left *= right;
    return left;
}

/*!
 * \fn inline BlockDiagonalMatrix operator*(double scaling, const BlockDiagonalMatrix & right)
 * \brief Multiplication operator: uniform scaling of BlockDiagonalMatrix version
 * \param scaling the scaling factor
 * \param object the square block diagonal matrix to be scaled
 */
template <typename T, int nrBlocks_, int dimBlock_>
inline BlockDiagonalMatrix<T, nrBlocks_, dimBlock_> operator*(double scaling,
        BlockDiagonalMatrix<T, nrBlocks_, dimBlock_> & object)
{
    object *= scaling;
    return object;
}

#endif // BLOCKDIAGONALMATRIX_HPP
