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

    module pedra_dblas

    use pedra_precision

    implicit none

    public dasum
    public daxpy
    public dscal
    public dswap
    public dcopy
    public dzero
    public dnorm2
    public idamax
    public vector_product

    contains

    real(kind=dp) function dasum(n, dx, incx)
!
! - Reference BLAS level1 routine (version 3.4.0) --
! - Reference BLAS is a software package provided by Univ. of Tennessee,    --
! - Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!   November 2011
!
!   .. Scalar Arguments ..
    integer(kind=regint_k), intent(in) :: incx, n
!   ..
!   .. Array Arguments ..
    real(kind=dp) :: dx(*)
!   ..
!
! ====================================================================
!
!   .. Local Scalars ..
    real(kind=dp) :: dtemp
    integer(kind=regint_k) :: i, m, mp1, nincx
!   ..
!   .. Intrinsic Functions ..
    intrinsic dabs, mod
!   ..
    dasum = 0.0d0
    dtemp = 0.0d0
    if (n.le.0 .or. incx.le.0) return
    if (incx.eq.1) then
!      code for increment equal to 1
!
!
!      clean-up loop
!
       m = mod(n, 6_regint_k)
       if (m.ne.0) then
          do i = 1,m
             dtemp = dtemp + dabs(dx(i))
          end do
          if (n.lt.6) then
             dasum = dtemp
             return
          end if
       end if
       mp1 = m + 1
       do i = mp1,n,6
          dtemp = dtemp + dabs(dx(i)) + dabs(dx(i+1)) +    &
                 dabs(dx(i+2)) + dabs(dx(i+3)) +           &
                 dabs(dx(i+4)) + dabs(dx(i+5))
       end do
    else
!
!      code for increment not equal to 1
!
       nincx = n*incx
       do i = 1,nincx,incx
          dtemp = dtemp + dabs(dx(i))
       end do
    end if
    dasum = dtemp

    end function dasum

    subroutine daxpy(n, da, dx, incx, dy, incy)
!
! - Reference BLAS level1 routine (version 3.4.0) --
! - Reference BLAS is a software package provided by Univ. of Tennessee,    --
! - Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!   November 2011
!
!   .. Scalar Arguments ..
    real(kind=dp), intent(in) :: da
    integer(kind=regint_k), intent(in) :: incx, incy, n
!   ..
!   .. Array Arguments ..
    real(kind=dp) :: dx(*), dy(*)
!   ..
!
! ====================================================================
!
!   .. Local Scalars ..
    integer(kind=regint_k) :: i, ix, iy, m, mp1
!   ..
!   .. Intrinsic Functions ..
    intrinsic mod
!   ..
    if (n.le.0) return
    if (da.eq.0.0D0) return
    if (incx.eq.1 .and. incy.eq.1) then
!
!      code for both increments equal to 1
!
!
!      clean-up loop
!
       m = mod(n, 4_regint_k)
       if (m.ne.0) then
          do i = 1,m
             dy(i) = dy(i) + da*dx(i)
          end do
       end if
       if (n.lt.4) return
       mp1 = m + 1
       do i = mp1,n,4
          dy(i) = dy(i) + da*dx(i)
          dy(i+1) = dy(i+1) + da*dx(i+1)
          dy(i+2) = dy(i+2) + da*dx(i+2)
          dy(i+3) = dy(i+3) + da*dx(i+3)
       end do
    else
!
!      code for unequal increments or equal increments
!        not equal to 1
!
       ix = 1
       iy = 1
       if (incx.lt.0) ix = (-n+1)*incx + 1
       if (incy.lt.0) iy = (-n+1)*incy + 1
       do i = 1,n
        dy(iy) = dy(iy) + da*dx(ix)
        ix = ix + incx
        iy = iy + incy
       end do
    end if

    end subroutine daxpy

    subroutine dzero(dx, length)

!...................................................................
! Last revision 5-May-1984 by Hans Jorgen Aa. Jensen
!
!   Subroutine DZERO sets a real array of length *LENGTH*
!   to zero.
!...................................................................
    integer(kind=regint_k), intent(in)    :: length
    real(kind=dp), intent(inout) :: dx(length)

    integer(kind=regint_k) :: i

    if (length <= 0) return

    do i = 1, length
        dx(i) = 0.0d0
    end do

    end subroutine dzero

    real(kind=dp) function dnorm2(n, x, incx)
!
! - Reference BLAS level1 routine (version 3.4.0) --
! - Reference BLAS is a software package provided by Univ. of Tennessee,    --
! - Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!   November 2011
!
!   .. Scalar Arguments ..
    integer(kind=regint_k), intent(in) :: incx, n
!   ..
!   .. Array Arguments ..
    real(kind=dp) :: x(*)
!   ..
!
! ====================================================================
!
!   .. Parameters ..
    real(kind=dp) :: one, zero
    parameter (one=1.0d+0,zero=0.0d+0)
!   ..
!   .. Local Scalars ..
    real(kind=dp) :: absxi, norm, scale, ssq
    integer(kind=regint_k) :: ix
!   ..
!   .. Intrinsic Functions ..
    intrinsic abs, sqrt
!   ..
    if (n.lt.1 .or. incx.lt.1) then
        norm = zero
    else if (n.eq.1) then
        norm = abs(x(1))
    else
        scale = zero
        ssq = one
!      The following loop is equivalent to this call to the LAPACK
!      auxiliary routine:
!      CALL DLASSQ( N, X, INCX, SCALE, SSQ )
!
        do 10 ix = 1,1 + (n-1)*incx,incx
            if (x(ix).ne.zero) then
                absxi = abs(x(ix))
                if (scale.lt.absxi) then
                    ssq = one + ssq* (scale/absxi)**2
                    scale = absxi
                else
                    ssq = ssq + (absxi/scale)**2
                end if
            end if
  10     continue
        norm = scale*sqrt(ssq)
    end if
!
    dnorm2 = norm

    end function dnorm2

    subroutine dscal(n, da, dx, incx)
!
! - Reference BLAS level1 routine (version 3.4.0) --
! - Reference BLAS is a software package provided by Univ. of Tennessee,    --
! - Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!   November 2011
!
!   .. Scalar Arguments ..
    real(kind=dp), intent(in) :: da
    integer(kind=regint_k), intent(in) :: incx, n
!   ..
!   .. Array Arguments ..
    real(kind=dp), intent(inout) :: dx(*)
!   ..
!
! ====================================================================
!
!   .. Local Scalars ..
    integer(kind=regint_k) :: i, m, mp1, nincx
!   ..
!   .. Intrinsic Functions ..
    intrinsic mod
!   ..
    if (n.le.0 .or. incx.le.0) return
    if (incx.eq.1) then
!
!      code for increment equal to 1
!
!
!      clean-up loop
!
       m = mod(n, 5_regint_k)
       if (m.ne.0) then
          do i = 1,m
             dx(i) = da*dx(i)
          end do
          if (n.lt.5) return
       end if
       mp1 = m + 1
       do i = mp1,n,5
          dx(i) = da*dx(i)
          dx(i+1) = da*dx(i+1)
          dx(i+2) = da*dx(i+2)
          dx(i+3) = da*dx(i+3)
          dx(i+4) = da*dx(i+4)
       end do
    else
!
!      code for increment not equal to 1
!
       nincx = n*incx
       do i = 1,nincx,incx
          dx(i) = da*dx(i)
       end do
    end if

    end subroutine dscal

    integer(kind=regint_k) function idamax(n, dx, incx)
!
! - Reference BLAS level1 routine (version 3.4.0) --
! - Reference BLAS is a software package provided by Univ. of Tennessee,    --
! - Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!   November 2011
!
!   .. Scalar Arguments ..
    integer(kind=regint_k), intent(in) :: incx, n
!   ..
!   .. Array Arguments ..
    real(kind=dp), intent(in) :: dx(*)
!   ..
!
! ====================================================================
!
!   .. Local Scalars ..
    real(kind=dp) :: dmax
    integer(kind=regint_k) :: i, ix
!   ..
!   .. Intrinsic Functions ..
    intrinsic dabs
!   ..
    idamax = 0
    if (n.lt.1 .or. incx.le.0) return
    idamax = 1
    if (n.eq.1) return
    if (incx.eq.1) then
!
!      code for increment equal to 1
!
       dmax = dabs(dx(1))
       do i = 2,n
          if (dabs(dx(i)).gt.dmax) then
             idamax = i
             dmax = dabs(dx(i))
          end if
       end do
    else
!
!      code for increment not equal to 1
!
       ix = 1
       dmax = dabs(dx(1))
       ix = ix + incx
       do i = 2,n
          if (dabs(dx(ix)).gt.dmax) then
             idamax = i
             dmax = dabs(dx(ix))
          end if
          ix = ix + incx
       end do
    end if

    end function idamax

    subroutine dswap(n, dx, incx, dy, incy)
!
!  -Reference BLAS level1 routine (version 3.4.0) --
!  -Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!   November 2011
!
!   .. Scalar Arguments ..
    integer(kind=regint_k), intent(in) :: incx, incy, n
!   ..
!   .. Array Arguments ..
    real(kind=dp) :: dx(*), dy(*)
!   ..
!
!  ===================================================================
!
!   .. Local Scalars ..
    real(kind=dp) :: dtemp
    integer(kind=regint_k) :: i, ix, iy, m, mp1
!   ..
!   .. Intrinsic Functions ..
    intrinsic mod
!   ..
    if (n.le.0) return
    if (incx.eq.1 .and. incy.eq.1) then
!
!     code for both increments equal to 1
!
!
!     clean-up loop
!
       m = mod(n, 3_regint_k)
       if (m.ne.0) then
          do i = 1,m
             dtemp = dx(i)
             dx(i) = dy(i)
             dy(i) = dtemp
          end do
          if (n.lt.3) return
       end if
       mp1 = m + 1
       do i = mp1,n,3
          dtemp = dx(i)
          dx(i) = dy(i)
          dy(i) = dtemp
          dtemp = dx(i+1)
          dx(i+1) = dy(i+1)
          dy(i+1) = dtemp
          dtemp = dx(i+2)
          dx(i+2) = dy(i+2)
          dy(i+2) = dtemp
       end do
    else
!
!     code for unequal increments or equal increments not equal
!       to 1
!
       ix = 1
       iy = 1
       if (incx.lt.0) ix = (-n+1)*incx + 1
       if (incy.lt.0) iy = (-n+1)*incy + 1
       do i = 1,n
          dtemp = dx(ix)
          dx(ix) = dy(iy)
          dy(iy) = dtemp
          ix = ix + incx
          iy = iy + incy
       end do
    end if

    end subroutine dswap

    subroutine dcopy(n, dx, incx, dy, incy)
!
! - Reference BLAS level1 routine (version 3.4.0) --
! - Reference BLAS is a software package provided by Univ. of Tennessee,    --
! - Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!   November 2011
!
!   .. Scalar Arguments ..
    integer(kind=regint_k), intent(in) :: incx, incy, n
!   ..
!   .. Array Arguments ..
    real(kind=dp) :: dx(*), dy(*)
!   ..
!
! ====================================================================
!
!   .. Local Scalars ..
    integer(kind=regint_k) :: i, ix, iy, m, mp1
!   ..
!   .. Intrinsic Functions ..
    intrinsic mod
!   ..
    if (n.le.0) return
    if (incx.eq.1 .and. incy.eq.1) then
!
!      code for both increments equal to 1
!
!
!      clean-up loop
!
       m = mod(n, 7_regint_k)
       if (m.ne.0) then
          do i = 1,m
             dy(i) = dx(i)
          end do
          if (n.lt.7) return
       end if
       mp1 = m + 1
       do i = mp1,n,7
          dy(i) = dx(i)
          dy(i+1) = dx(i+1)
          dy(i+2) = dx(i+2)
          dy(i+3) = dx(i+3)
          dy(i+4) = dx(i+4)
          dy(i+5) = dx(i+5)
          dy(i+6) = dx(i+6)
       end do
    else
!
!      code for unequal increments or equal increments
!        not equal to 1
!
       ix = 1
       iy = 1
       if (incx.lt.0) ix = (-n+1)*incx + 1
       if (incy.lt.0) iy = (-n+1)*incy + 1
       do i = 1,n
          dy(iy) = dx(ix)
          ix = ix + incx
          iy = iy + incy
       end do
    end if

    end subroutine dcopy

    subroutine vector_product(p1, p2, p3, dnorm3)
! Calculates vector product and norm of resulting vector

    real(kind=dp), intent(in)  :: p1(3), p2(3)
    real(kind=dp), intent(out) ::p3(3), dnorm3

    p3(1) = p1(2)*p2(3) - p1(3)*p2(2)
    p3(2) = p1(3)*p2(1) - p1(1)*p2(3)
    p3(3) = p1(1)*p2(2) - p1(2)*p2(1)
    dnorm3 = sqrt(p3(1)*p3(1) + p3(2)*p3(2) + p3(3)*p3(3))

    end subroutine vector_product

    end module pedra_dblas
