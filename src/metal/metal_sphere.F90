!pcmsolver_copyright_start
!      PCMSolver, an API for the Polarizable Continuum Model
!      Copyright (C) 2013 Roberto Di Remigio, Luca Frediani and contributors
!      
!      This file is part of PCMSolver.
!
!      PCMSolver is free software: you can redistribute it and/or modify       
!      it under the terms of the GNU Lesser General Public License as published by
!      the Free Software Foundation, either version 3 of the License, or
!      (at your option) any later version.
!                                                                           
!      PCMSolver is distributed in the hope that it will be useful,
!      but WITHOUT ANY WARRANTY; without even the implied warranty of
!      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!      GNU Lesser General Public License for more details.
!                                                                           
!      You should have received a copy of the GNU Lesser General Public License
!      along with PCMSolver.  If not, see <http://www.gnu.org/licenses/>.
!
!      For information on the complete list of contributors to the
!      PCMSolver API, see: <https://repo.ctcc.no/projects/pcmsolver>
!pcmsolver_copyright_end

    module metal_sphere

    use, intrinsic :: iso_c_binding

    implicit none

    public greens_function 
    
    private

    contains
      
    subroutine greens_function(epssol, epsre, epsim, radius,   &
         ps, p1, p2, greenre, greenim) bind(c, name='greens_function')

! Passed variables      
    real(c_double), intent(in)  :: epssol, epsre, epsim, radius
    real(c_double), intent(in)  :: p1(*), p2(*), ps(*)
    real(c_double), intent(out) :: greenre, greenim
! Local variables      
    complex(16) :: green, eps2, ui
      
    ui = (0.0, 1.0)                                               
    eps2 = epsre + epsim * ui                                     
    green = gsfera(epssol, eps2, ps(1), ps(2), ps(3), radius,  &  
                  p1(1), p1(2), p1(3),                         &  
                  p2(1), p2(2), p2(3))                            
    greenre = real(real(green))
    greenim = real(aimag(green))

    end subroutine greens_function                                

    complex(16) function gsfera(eps, eps2, xs, ys, zs, rs,     &
                                xi, yi, zi, xj, yj, zj)
! Passed variables 
    real(8), intent(in)     :: eps
    complex(16), intent(in) :: eps2
    real(8), intent(in)     :: xs, ys, zs, rs ! Sphere center and radius 
    real(8), intent(in)     :: xi, yi, zi     ! Source point
    real(8), intent(in)     :: xj, yj, zj     ! Probe point
! Local variables          
    complex(16) :: coefl
    real(8)     :: di, dj
    real(8)     :: dim, xim, yim, zim, qim    ! Image quantities
    real(8)     :: gc, cost, aa, arg, argl, gsfera2
    integer     :: maxl, l, m

!lf print *, eps,eps2
!lf print *, 'Sphere', xs,ys,zs,rs
!lf print *, 'Source', xi,yi,zi
!lf print *, 'Probe ', xj,yj,zj
    di = sqrt((xi - xs)**2 + (yi - ys)**2 + (zi - zs)**2)
    dj = sqrt((xj - xs)**2 + (yj - ys)**2 + (zj - zs)**2)
    dim = rs**2 / di
    xim = dim * (xi - xs) / (di + xs) 
    yim = dim * (yi - ys) / (di + ys) 
    zim = dim * (zi - zs) / (di + zs) 
    qim = rs / di
    gc = qim / sqrt((xj - xim)**2 + (yj - yim)**2 + (zj - zim)**2)
    gc = gc - qim / sqrt((xj - xs)**2 + (yj - ys)**2 + (zj - zs)**2)
    cost = (xj - xs) * (xi - xs) + (yj - ys) * (yi - ys) + (zj - zs) * (zi - zs)
    cost = cost / (di * dj)
    aa = (rs * rs / (di * dj))**4.
    maxl = int(800. / (abs(eps2 + eps))**0.4 * aa)
    maxl = max(maxl, 10)
!lf write (6,*) "maxl",maxl
    arg = rs * rs / (di * dj)
    argl = rs / (di * dj)
    gsfera2 = 0.0d0
    m = 0
    do l = 1, maxl
      argl = argl * arg
      coefl = (eps2 - eps) * float(l) / ((eps2 + eps) * float(l) + 1.0d0)
      coefl = coefl - (eps2 - eps) / (eps2 + eps)
      gsfera = gsfera + coefl * argl * legendre_polynomial(l, m, cost)
    enddo
    gsfera = gsfera + gc * (eps2- eps) / (eps2 + eps)
!lf: WARNING: sign change and divided by epsilon of solvent!!!!!!
    gsfera = -gsfera / eps
! DEBUGGING STUFF      
!   arg = rs * rs / (di * dj)
!   argl = rs / (di * dj)
!   gsfera = 0.0d0
!   m = 0
!   do l = 1, 5000
!     argl = argl * arg
!     coefl = (eps2 - eps) * float(l) / ((eps2 + eps) * float(l) + 1.0d0)
!     gsfera = gsfera + coefl * argl * legendre_polynomial(l, m, cost)
!   enddo 
!lf write (6, *) "gsfera2", gsfera
    end function gsfera
    
    real(8) function legendre_polynomial(l, m, x)
! Computes the associated Legendre polynomial P_l^m(x)
! Passed variables      
    integer, intent(in) :: l, m
    real(8), intent(in) :: x
! Local variables      
    integer :: i, ll
    real(8) :: fact, pll, pmm, pmmp1, somx2

    pmm = 1.0d0
    if (m .gt. 0) then
      somx2 = sqrt((1.0 - x) * (1.0 + x))
      fact = 1.0d0
      do i = 1, m
        pmm = -pmm * fact * somx2
        fact = fact + 2.0d0
      enddo
    endif
    if (l .eq. m) then
      legendre_polynomial = pmm
    else
      pmmp1 = x * (2 * m + 1) * pmm
      if (l .eq. m + 1) then
        legendre_polynomial = pmmp1
      else
        do ll = m + 2, l
          pll = (x * (2.0 * ll - 1) * pmmp1 - (ll + m - 1) * pmm) / (ll - m)
          pmm = pmmp1
          pmmp1 = pll
        enddo
        legendre_polynomial = pll
      endif
    endif

    end function legendre_polynomial

    end module metal_sphere 
