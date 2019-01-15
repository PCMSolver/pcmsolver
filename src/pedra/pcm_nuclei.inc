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

      integer(kind=regint_k) NUCPRE, NUCNUM, NUCDEG, ISTBNU, NCTOT, &
             NUCIND, NUCDEP, NTRACO, ITRACO, NATOMS, NFLOAT,        &
             NBASIS, NLARGE, NSMALL, NPBAS,  NPLRG,  NPSML,         &
             NCHTOT, INCENT, INUNIQ, NDEGNM, ISOTOP, IZATOM,        &
             NBASISAUX, NPBASAUX, NAUX, NPAUX
#if defined (SYS_CRAY)
      REAL    CHARGE, CORD, GNUEXP
#else
      DOUBLE PRECISION    CHARGE, CORD, GNUEXP
#endif
      LOGICAL NOORBT,GAUNUC
      COMMON /PCM_NUCLEI/ CHARGE(MXCENT), CORD(3,MXCENT),GNUEXP(MXCENT),       &
                     GAUNUC, NOORBT(MXCENT), NUCPRE(MXCENT),                   &
                     NUCNUM(MXCENT,8), NUCDEG(MXCENT), ISTBNU(MXCENT),         &
                     NDEGNM(MXCENT), NUCIND, NUCDEP, NTRACO, ITRACO(3),        &
                     NATOMS, NFLOAT, NBASIS, NLARGE, NSMALL, NPBAS,            &
                     NPLRG, NPSML, NCHTOT, INCENT(MXCENT), NCTOT,              &
                     INUNIQ(MXCENT), ISOTOP(MXCENT),IZATOM(MXCENT),            &
                     NBASISAUX, NPBASAUX, NAUX, NPAUX
      CHARACTER NAMEX*6, NAMDEP*6, NAMDPX*8, NAMN*4
      COMMON /PCM_NUCLEC/ NAMEX(MXCOOR), NAMDEP(MXCENT), NAMDPX(MXCOOR),       &
                     NAMN(MXCENT)
      INTEGER(kind=regint_k) MULBSI
      COMMON /PCM_MULBAS/ MULBSI(MXCENT)
!     MULBAS has been added for multiple basis sets (WK/UniKA/31-10-2002).
