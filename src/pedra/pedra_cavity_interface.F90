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

!
!     simple input reader for cavity generator
!     written by Krzysztof Mozgawa, 2010
!
!     RDR, 280114. Put things in makecav.F inside here directly.
!
subroutine generatecavity_cpp(maxts_, maxsph_, maxvert_,         &
    xtscor_, ytscor_, ztscor_, ar_,                              &
    xsphcor_, ysphcor_, zsphcor_, rsph_,                         &
    nts_, ntsirr_, nesfp_, addsph_,                              &
    xe_, ye_, ze_, rin_, masses_, avgArea_, rsolv_, ret_,        &
    nr_gen_, gen1_, gen2_, gen3_,                                &
    nvert_, vert_, centr_, isphe_, pedra_, len_pedra_)           &
    bind(C, name='generatecavity_cpp')

use, intrinsic :: iso_c_binding
use pedra_precision
use pedra_symmetry, only: get_point_group, point_group
use pedra_cavity, only: polyhedra_driver
use strings, only: carray_to_fstring

implicit none

#include "pcm_pcmdef.inc"
#include "pcm_mxcent.inc"
#include "pcm_pcm.inc"

integer(c_int)  :: maxts_, maxsph_, maxvert_
real(c_double)  :: xtscor_(maxts_), ytscor_(maxts_), ztscor_(maxts_)
real(c_double)  :: xsphcor_(maxts_), ysphcor_(maxts_), zsphcor_(maxts_), rsph_(maxts_)
real(c_double)  :: ar_(maxts_), xe_(maxts_), ye_(maxts_), ze_(maxts_), rin_(maxts_)
real(c_double)  :: masses_(maxts_)
real(c_double)  :: avgArea_, rsolv_, ret_
integer(c_int)  :: nts_, ntsirr_, nesfp_, addsph_
integer(c_int)  :: nr_gen_, gen1_, gen2_, gen3_
integer(c_int)  :: nvert_(maxts_)
real(c_double)  :: vert_(maxts_ * 30), centr_(maxts_ * 30)
integer(c_int)  :: isphe_(maxts_)
integer(c_int) :: len_pedra_
character(kind=c_char, len=1), intent(in) :: pedra_(len_pedra_+1)

integer(c_int)    :: i, j, k, offset
integer(c_int)    :: error_code
integer(kind=regint_k) :: pedra_unit
logical           :: pedra_open, pedra_exist
real(c_double), allocatable :: vert(:, :, :), centr(:, :, :)
character(len=len_pedra_) :: pedra
type(point_group) :: pgroup


pedra = carray_to_fstring(pedra_)
!lvpri = 121201_regint_k
! The following INQUIRE statement returns whether the file named cavity.off is
! connected in logical variable off_open, whether the file exists in logical
! variable off_exist, and the unit number in integer(kind=regint_k) variable off_unit
pedra_unit = 121201_regint_k
inquire(file = pedra, opened = pedra_open, &
        exist = pedra_exist)
if (pedra_exist) then
   open(unit = pedra_unit,                   &
       file = pedra,       &
       status = 'unknown',       &
       form = 'formatted',       &
       access = 'sequential')
   close(unit = pedra_unit, status = 'delete')
end if
open(unit = pedra_unit,                      &
    file = pedra,          &
    status = 'new',              &
    form = 'formatted',          &
    access = 'sequential')
rewind(pedra_unit)
areats = avgArea_
icesph = 1
iprpcm = 3
call get_point_group(pedra_unit, pgroup, nr_gen_, gen1_, gen2_, gen3_)
rsolv = rsolv_
! These parameters are fixed see one of the original GePol papers
omega = 40.0d+00
fro = 0.7d+00
! ret is the minimum radius of added spheres
ret = ret_
! nesfp is the number of initial spheres.
! nesf is the total number of spheres: initial + added
nesfp = nesfp_
do i = 1, nesfp
   xe(i) = xe_(i)
   ye(i) = ye_(i)
   ze(i) = ze_(i)
   rin(i) = rin_(i)
   alpha(i) = 1.0d0
end do
! Workaround to silence compilation warning
maxsph_ = mxsp
maxvert_= mxver

! Allocate space for the arrays containing the vertices and the centers
! of the tesserae arcs
allocate(vert(mxts, 10, 3))
vert = 0.0d0
allocate(centr(mxts, 10, 3))
centr = 0.0d0

nesf = nesfp

call polyhedra_driver(pgroup, vert, centr, masses_, pedra_unit, error_code)

! Common block dark magic, it will disappear one day...
nts_ = nts
ntsirr_ = ntsirr
! Pass the number of added spheres back, to update the GePolCavity
! object in the right way.
! nesf: total number of spheres; nesfp: number of original spheres
addsph_ = nesf - nesfp
do i = 1, nts
   xtscor_(i) = xtscor(i)
   ytscor_(i) = ytscor(i)
   ztscor_(i) = ztscor(i)
   ar_(i) = as(i)
   xsphcor_(i) = xe(isphe(i))
   ysphcor_(i) = ye(isphe(i))
   zsphcor_(i) = ze(isphe(i))
   rsph_(i) = re(isphe(i))
   nvert_(i) = nvert(i)
! Fill the vert_ and centr_ arrays
   do j = 1, nvert(i)
      do k = 1, 3
          offset = i + j * nts + k * nts * nvert(i)
          vert_(offset) = vert(i, j, k)
          centr_(offset) = centr(i, j, k)
      end do
   end do
end do

do i = 1, nesf
  isphe_(i) = isphe(i)
  xe_(i)  = xe(i)
  ye_(i)  = ye(i)
  ze_(i)  = ze(i)
  rin_(i) = re(i)
end do

! Clean-up
deallocate(vert)
deallocate(centr)

write(pedra_unit, *) "Error code is ", error_code
write(pedra_unit, *) '<<< Done with PEDRA Fortran code >>>'

close(pedra_unit)

end subroutine generatecavity_cpp
