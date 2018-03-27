!
! PCMSolver, an API for the Polarizable Continuum Model
! Copyright (C) 2018 Roberto Di Remigio, Luca Frediani and contributors.
!
! This file is part of PCMSolver.
!
! PCMSolver is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! PCMSolver is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with PCMSolver.  If not, see <http://www.gnu.org/licenses/>.
!
! For information on the complete list of contributors to the
! PCMSolver API, see: <http://pcmsolver.readthedocs.io/>
!

module tsless_lapack

use tsless_precision

implicit none

public solve_linear_system

private

contains

    !> \brief Solve linear system of equations
    !> \author Roberto Di Remigio
    !> \year 2015
    !> \param[in] ndim dimension of the problem
    !> \param[in] A system matrix
    !> \param[in] b right-hand side
    !> \param[out] x solution vector
    !> Computes the explicit inverse by LU decomposition and then
    !> performs matrix-vector multiplication A^-1 b = x
    pure subroutine solve_linear_system(ndim, A, b, x)

    !> Passed variables
    integer(kind=regint_k), intent(in) :: ndim
    real(kind=dp),       intent(inout) :: A(ndim, ndim)
    real(kind=dp),          intent(in) :: b(ndim)
    real(kind=dp),         intent(out) :: x(ndim)
    !> Local variables
    real(kind=dp), allocatable :: Ainverse(:, :)
    integer(kind=regint_k) :: i, j

    !> Compute inverse by LU decomposition
    allocate(Ainverse(ndim, ndim)); Ainverse = 0.0_dp
    call inverse_LU(ndim, A, Ainverse)

    !> Solve linear system
    x = 0.0_dp
    do i = 1_regint_k, ndim
       do j = 1_regint_k, ndim
           x(i) = x(i) + Ainverse(i, j) * b(j)
       end do
    end do

    deallocate(Ainverse)

    end subroutine solve_linear_system

    !> \brief Computes inverse matrix
    !> \author Alexander L. Godunov
    !> \year 2009
    !> \param[in, out] A matrix to be inverted
    !> \param[in, out] C inverse matrix
    !> \param[in] ndim dimension of the matrix
    !> The method is based on Doolittle LU factorization for Ax = b
    !> notice that the contents of the original matrix will be destroyed
    pure subroutine inverse_LU(ndim, A, C)

    !> Passed variables
    integer(kind=regint_k), intent(in) :: ndim
    real(kind=dp),       intent(inout) :: A(ndim, ndim)
    real(kind=dp),       intent(inout) :: C(ndim, ndim)
    !> Local variables
    real(kind=dp), allocatable :: L(:, :), U(:, :)
    real(kind=dp), allocatable :: b(:), d(:), x(:)
    real(kind=dp) :: coeff
    integer(kind=regint_k) :: i, j, k

    ! step 0: initialization for matrices L and U and b
    allocate(L(ndim, ndim)); L = 0.0_dp
    allocate(U(ndim, ndim)); U = 0.0_dp
    allocate(b(ndim)); b = 0.0_dp
    allocate(d(ndim)); d = 0.0_dp
    allocate(x(ndim)); x = 0.0_dp

    ! step 1: forward elimination
    do k = 1, ndim-1
        do i = k + 1, ndim
            coeff = A(i, k) / A(k, k)
            L(i, k) = coeff
            do j = k + 1, ndim
                A(i, j) = A(i, j) - coeff * A(k, j)
            end do
        end do
    end do

    ! Step 2: prepare L and U matrices
    ! L matrix is a matrix of the elimination coefficient
    ! + the diagonal elements are 1.0
    do i = 1, ndim
        L(i, i) = 1.0_dp
    end do
    ! U matrix is the upper triangular part of A
    do j = 1, ndim
        do i = 1, j
            U(i, j) = A(i, j)
        end do
    end do

    ! Step 3: compute columns of the inverse matrix C
    do k = 1, ndim
        b(k) = 1.0_dp
        d(1) = b(1)
        ! Step 3a: Solve Ld=b using the forward substitution
        do i = 2, ndim
            d(i) = b(i)
            do j = 1, i - 1
                d(i) = d(i) - L(i, j) * d(j)
            end do
        end do
        ! Step 3b: Solve Ux=d using the back substitution
        x(ndim) = d(ndim) / U(ndim, ndim)
        do i = ndim - 1, 1, -1
            x(i) = d(i)
            do j = ndim, i + 1, -1
                x(i) = x(i) - U(i, j) * x(j)
            end do
            x(i) = x(i) / U(i, i)
        end do
        ! Step 3c: fill the solutions x(ndim) into column k of C
        do i = 1, ndim
            C(i, k) = x(i)
        end do
        b(k) = 0.0_dp
    end do

    deallocate(L)
    deallocate(U)
    deallocate(b)
    deallocate(d)
    deallocate(x)

    end subroutine inverse_LU

end module tsless_lapack
