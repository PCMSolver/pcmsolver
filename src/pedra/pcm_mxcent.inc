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
!     file: mxcent.h
!
!     MXCENT = max number of nuclei + point charges + ghost orbital centers
!
!     IF you change MXCENT you should do a "make depend"
!     and then rebuild the program using the command "make".
!
!     In case of QM3 MXNEW is used to allocate memory in herrdn.F.
!     To run a QM3 calculation in most cases MXCENT will
!     have to be around 2000 - 3000!!! Remember to set MXQM3 = MXCENT in
!     qm3.h!!!
!
      integer(kind=regint_k) MXNEW, MXCENT, MXCOOR
      PARAMETER (MXNEW =1000, MXCENT = 1000, MXCOOR = 3*MXCENT)
