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

module pedra_precision
! Read this: http://stackoverflow.com/a/3170438/2528668
! and this:  http://stackoverflow.com/a/3204981/2528668

implicit none

! Integer types
! 32-bit integers
integer, parameter :: regint_k   = selected_int_kind(8)
! 64-bit integers
integer, parameter :: largeint_k = selected_int_kind(18)

! Real types
! Single-precision real
integer, parameter :: sp = kind(1.0)
! Double-precision real
integer, parameter :: dp = selected_real_kind(2*precision(1.0_sp))

end module pedra_precision
