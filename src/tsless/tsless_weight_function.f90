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

module tsless_weight_function

use, intrinsic :: iso_c_binding
use tsless_precision

implicit none

!> \brief Describes a weight function
!> Data fiels have dimensions 2*nr_derivative+3
type, public :: weight_function
    !> Number of continuous derivatives
    integer(kind=regint_k) :: nr_derivative
    !>
    real(kind=dp), pointer :: weight_1(:)
    !>
    real(kind=dp), pointer :: weight_2(:)
    !>
    real(kind=dp), pointer :: weight_3(:)
end type weight_function

public create_weight_function
public destroy_weight_function
public scaling_factor

private

contains

    function create_weight_function(nderiv, normalization, dmin, printer) result(wfun)

    use, intrinsic :: iso_fortran_env, only: output_unit

    !> Passed variables
    integer(kind=regint_k), intent(in) :: nderiv
    integer(kind=regint_k), intent(in) :: normalization
    real(kind=dp),          intent(in) :: dmin
    integer, optional,      intent(in) :: printer
    !> Output variables
    type(weight_function) :: wfun
    !> Local variables
    integer :: print_out

    if (present(printer)) then
       print_out = printer
    else
       print_out = output_unit
    end if

    call allocate_weight_function(wfun, nderiv)
    call compute(wfun, normalization, dmin, print_out)

    end function create_weight_function

    pure subroutine allocate_weight_function(wfun, nderiv)

    !> Passed variables
    type(weight_function), intent(inout) :: wfun
    integer(kind=regint_k), intent(in)   :: nderiv
    !> Local variables
    integer(kind=regint_k) :: nparms

    wfun%nr_derivative = nderiv
    nparms = 2_regint_k * nderiv + 3_regint_k
    allocate(wfun%weight_1(nparms)); wfun%weight_1 = 0.0_dp
    allocate(wfun%weight_2(nparms)); wfun%weight_2 = 0.0_dp
    allocate(wfun%weight_3(nparms)); wfun%weight_3 = 0.0_dp

    end subroutine allocate_weight_function

    pure subroutine destroy_weight_function(wfun)

    !> Passed variables
    type(weight_function), intent(inout) :: wfun

    wfun%nr_derivative = 0_regint_k
    deallocate(wfun%weight_1)
    deallocate(wfun%weight_2)
    deallocate(wfun%weight_3)

    end subroutine destroy_weight_function

    !> \brief
    !> \author Christian S. Pomelli
    !> \param[inout] wfun weight function
    !> \param[in] nord maximum number of continuous derivatives
    !> \param[in] ifun
    !> \param[in] dmin minimal distance between sampling points
    !> \param[in] printer logical unit for printing
    !> ifun = 0 non-normalized function (gamma)
    !> ifun = 1 normalized function (omega)
    subroutine compute(wfun, ifun, dmin, printer)

    use, intrinsic :: iso_fortran_env, only: output_unit

    use tsless_lapack, only: solve_linear_system

    !> Passed variables
    type(weight_function), intent(inout) :: wfun
    integer(kind=regint_k),   intent(in) :: ifun
    real(kind=dp),            intent(in) :: dmin
    integer, optional,        intent(in) :: printer
    !> Local variables
    real(kind=dp), allocatable :: C(:, :)
    real(kind=dp), allocatable :: b(:)
    integer :: i, j, k, n, m
    integer :: print_out

    if (present(printer)) then
       print_out = printer
    else
       print_out = output_unit
    end if

    ! check parameters
    if ((ifun .lt. 0) .or. (ifun .gt. 1)) Then
       write(print_out, *) 'Unknown IFUN'
       stop
    end if

    ! m = number of parameters
    n = wfun%nr_derivative
    if (ifun .eq. 0) then
       m = 2*n + 2
    else
       m = 2*n + 3
    end if
    !> Allocate, clean up and fill coefficient matrix C
    allocate(C(m, m)); C = 0.0_dp
    !   l = 1 boundary conditions
    do j = 0_regint_k, n
       do k = j, m - 1_regint_k
          C(j+1, k+1) = coefficients(j, k)
       end do
    end do
    !   l = d boundary conditions
    do j = 0_regint_k, n
       do k = j, m - 1_regint_k
          C(j+2+n, k+1) = coefficients(j, k)
          C(j+2+n, k+1) = C(j+2+n, k+1) * dmin**(k-j)
       end do
    end do
    !   normalization condition
    if (ifun .eq. 1_regint_k) then
       do k = 0_regint_k, m - 1_regint_k
          C(m, k+1) = (1.0_dp - dmin**(k+1)) / (1.0_dp * (k+1))
       end do
    end if
    ! Allocate, clean up and fill b vector (RHS)
    allocate(b(m)); b = 0.0_dp
    b(1) = 1.0_dp
    if (ifun .eq. 1) b(m) = 1.0_dp - dmin
    ! 2. Solve linear system of equations
    call solve_linear_system(ndim=m, A=C, b=b, x=wfun%weight_1)
    ! write function and derivatives
    do i = 1_regint_k, m - 1_regint_k
        wfun%weight_2(i) = i * wfun%weight_1(i+1)
    end do

    do i = 1_regint_k, m - 2_regint_k
        wfun%weight_3(i) = i * wfun%weight_2(i+1)
    end do

    write(print_out, *) 'ifun = ', ifun
    do i = 1_regint_k, m
       write(print_out, '(a, i2, a, 3F15.4)') 'c(', i, ') = ', wfun%weight_1(i), wfun%weight_2(i), wfun%weight_3(i)
    end do
    write(print_out, *) '-------------------'

    deallocate(b)
    deallocate(C)

    end subroutine compute

    !> \brief Calculates scaling factor for weight of point
    !> \author Christian S. Pomelli
    !> \param[in] wfun weight function
    !> \param[in] arg argument of the regularized step function
    !> \param[in] dmin minimal distance between sampling points
    pure function scaling_factor(wfun, arg, dmin) result(scaling)

    !> Input variables
    type(weight_function), intent(in) :: wfun
    real(kind=dp),         intent(in) :: arg
    real(kind=dp),         intent(in) :: dmin
    !> Output variables
    real(kind=dp) :: scaling(3)
    !> Local variables
    integer :: i, nr_params
    real(kind=dp) :: tn

    ! Simply a polynomial. later add derivatives
    scaling = 0.0_dp
    tn = 1.0_dp
    ! number of polynomial parameters
    nr_params = 2_regint_k * wfun%nr_derivative + 3_regint_k
    do i = 1_regint_k, nr_params
        scaling(1) = scaling(1) + wfun%weight_1(i) * tn
        scaling(2) = scaling(2) + wfun%weight_2(i) * tn
        scaling(3) = scaling(3) + wfun%weight_3(i) * tn
        tn = tn * arg
    end do

    end function scaling_factor

    pure function coefficients(j, k) result(f)

    !> Passed variables
    integer(kind=regint_k), intent(in) :: j, k
    real(kind=dp) :: f
    !> Local variables
    integer :: i
    real(kind=dp) :: dkk

    f = 1.0_dp
    dkk = dble(k)
    if (j .ne. 0) then
       do i = 1, j
          f = f * dkk
          dkk = dkk - 1.0_dp
       end do
    end if

    end function coefficients

end module tsless_weight_function
