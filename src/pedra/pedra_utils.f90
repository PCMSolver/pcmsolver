!
! PCMSolver, an API for the Polarizable Continuum Model
! Copyright (C) 2019 Roberto Di Remigio, Luca Frediani and contributors.
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

module pedra_utils

  use pedra_precision
  use pedra_dlapack, only: dsyevj3, order

  implicit none

  public around
  public wlkdin

  private

contains

  subroutine around(head, print_unit)

  character(len=*), intent(in) :: head
  integer(kind=regint_k),          intent(in) :: print_unit

  integer(kind=regint_k) :: lhead, lng, ind, i

  lhead  = len_trim(head)
  lng    = lhead + 2
  ind = max(1,(80 - lng)/2 + 1)
  write (print_unit,'(//150a)') (' ',i=1,ind), '+', ('-',i=1,lng), '+'
  write (print_unit,'(150a)')   (' ',i=1,ind), '! ', head(1:lhead), ' !'
  write (print_unit,'(150a)')   (' ',i=1,ind), '+', ('-',i=1,lng), '+'
  write (print_unit,'()')

  end subroutine around

  ! WLKDIN on the basis of ABACUS abawalk.F
  ! written and/or modified by Krzysztof Mozgawa
  ! and Ville Weijo, 2010

  ! RDR 051213 revision. Remove call to LAPACK diagonalization
  ! subroutine. Replace with a custom function specialized
  ! for 3x3 matrices (thus avoiding having to link BLAS/LAPACK)
  subroutine wlkdin(cor, tmass, n, angmom, tinert, omega, cepval, cepvec, docopy, planar, linear)

  integer(kind=regint_k),    intent(in) :: n
  real(kind=dp),    intent(in) :: cor(n, 3)
  real(kind=dp),    intent(in) :: tmass(n)
  real(kind=dp),    intent(in) :: angmom(3)
  real(kind=dp), intent(inout) :: tinert(3, 3), cepvec(3, 3)
  real(kind=dp),   intent(out) :: omega(3), cepval(3)
  logical,   intent(out) :: planar, linear
  logical,    intent(in) :: docopy

  integer(kind=regint_k) :: i, j, k
  real(kind=dp) :: eigval(3), eigvec(3,3), tinver(3,3), eigvalinv(3,3), temp(3,3)
  ! Threshold to convert a numerical zero to a hard zero
  real(kind=dp), parameter :: tstlin = 1.0d-05
  real(kind=dp) :: average, r2

  tinert = 0.0d0
  eigvec = 0.0d0
  eigval = 0.0d0
  eigvalinv = 0.0d0
  tinver = 0.0d0
  temp = 0.0d0
  r2 = 0.0d0
  ! Build tensor of inertia:
  ! T_ij = sum_l m_l((x_l)_k(x_l)_k delta_ij - (x_l)_i(x_l)_j)
  do i = 1, n ! Loop on centers
     ! Calculate diagonal elements
     r2 = tmass(i) * dot_product(cor(i, :), cor(i, :))
     do j = 1, 3 ! Loop on coordinates
        ! Build diagonal elements
        tinert(j, j) = tinert(j, j) + r2
        do k = 1, 3 ! Loop on coordinates
           ! Build off-diagonal elements
           tinert(j, k) = tinert(j, k) - tmass(i)*cor(i,j)*cor(i,k)
        end do
     end do
  end do

  ! Symmetrize the inertia tensor
  average = 0.0d0
  do i = 1, 3
     do j = 1, 3
         average = 0.5d0 * (tinert(i, j) + tinert(j, i))
         tinert(i, j) = average
         tinert(j, i) = average
     end do
  end do

  temp = tinert
  ! Now diagonalize it, we want back both the eigenvalues and eigenvectors
  call dsyevj3(temp, eigvec, eigval)
  ! Eigenvalues and eigenvectors are sorted in decreasing order
  call order(eigvec, eigval, 3_regint_k, 3_regint_k)

  if ( abs(eigval(3)-eigval(2)-eigval(1)) .lt. tstlin) then
     planar = .true.
  else
     planar = .false.
  end if
  if (eigval(3) .lt. tstlin) then
     linear = .true.
     planar = .false.
     eigval(3) = 0.0d0
  else
     linear = .false.
     eigval(3) = 1d0/eigval(3)
  end if
  eigval(2) = 1.0d0/eigval(2)
  eigval(1) = 1.0d0/eigval(1)
  eigvalinv = 0.0d0
  eigvalinv(1,1) = eigval(1)
  eigvalinv(2,2) = eigval(2)
  eigvalinv(3,3) = eigval(3)

  tinver= matmul(matmul(eigvec, eigvalinv), transpose(eigvec))

  omega = matmul(tinver, angmom)

  if (docopy) then
     cepval = eigval
     cepvec = eigvec
  end if

  end subroutine wlkdin

end module pedra_utils

module strings

  use, intrinsic :: iso_c_binding, only: c_char, c_null_char

  implicit none

contains

  ! \brief Convert a Fortran string into a C string.
  ! \param[in] string_f03 a Fortran character string.
  ! \return array_c Null-terminated C string in a character array.
  pure function fstring_to_carray(string_f03) result(array_c)
    character(len=*), intent(in) :: string_f03
    character(kind=c_char, len=1) :: array_c(len(string_f03)+1)

    integer :: i

    do i = 1, len(string_f03)
        array_c(i) = string_f03(i:i)
    end do
    array_c(i) = c_null_char
  end function

  ! \brief Convert a C string into a Fortran string.
  ! \param[in] array_c a null-terminated C string in a character array.
  ! \return string_f03 Fortran character string (without null termination).
  pure function carray_to_fstring(array_c) result(string_f03)
    character(kind=c_char), intent(in) :: array_c(:)
    character(len=size(array_c)-1) :: string_f03

    integer :: i

    ! Don't copy any (presumably garbage) from beyond the null character.
    string_f03 = ''
    do i = 1, size(array_c)-1
        if (array_c(i) == c_null_char) exit
        string_f03(i:i) = array_c(i)
    end do
  end function

end module
