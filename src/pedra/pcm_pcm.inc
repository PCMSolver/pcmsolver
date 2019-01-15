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

! RET: minimum radius of added spheres
      DOUBLE PRECISION ALPHA, AREATS, AS, CCX, CCY, CCZ, DERCEN,   &
      DERRAD, DR, EPS, EPSINF, FN, FE, FRO, OMEGA, RET,        &
      RIN, RSOLV, QSN, QSE, QSENEQ, QLOC, PCMEN, PCMEE, PCMNE,     &
      PCMNN, QNUC, XE, YE, ZE, RE, STOT, VOL, RSPH,          &
      XTSCOR, YTSCOR, ZTSCOR, RDIF, POTCAVNUC, POTCAVELE
      integer(kind=regint_k) IAN, ICOMPCM, ICESPH, IDXSPH, INA, IPOLYG,           &
      IPRPCM, IRETCAV, ISPHE, NESF, NESFP, NC, NVERT, NTS, NTSIRR, &
      NRWCAV, LUPCMD, LUCAVD, LUPCMI, NPCMMT
      COMMON /PCM_CAV/ OMEGA,RET,FRO,ALPHA(MXSP),RIN(MXSP),ICESPH, &
                     IPRPCM,IRETCAV,IPOLYG,AREATS
      COMMON /PCM_DAT/ EPS,EPSINF,DR,RSOLV,ICOMPCM,NPCMMT
      COMMON /PCM_CH/  QSN(MXTS),QSE(MXTS),QSENEQ(MXTS),QLOC(MXTS,3), &
                     PCMEN,PCMEE,PCMNE,PCMNN,QNUC,FN,FE
      COMMON /PCM_PLY/ XE(MXSP),YE(MXSP),ZE(MXSP),RE(MXSP), &
                     ISPHE(MXTS),STOT,VOL,NESF,NESFP,NC(30)
      COMMON /PCM_SPH/ INA(MXSP),IAN(MXCENT), IDXSPH(MXCENT),          &
                      RSPH(MXCENT)
      COMMON /PCM_TES/ CCX,CCY,CCZ,XTSCOR(MXTSPT),YTSCOR(MXTSPT),      &
                     ZTSCOR(MXTSPT),AS(MXTS),RDIF,NVERT(MXTS),NTS,     &
                     NTSIRR,NRWCAV
      COMMON /PCM_LU/  LUPCMD,LUCAVD,LUPCMI(8)
      COMMON /PCM_DER/ DERCEN(MXSP,MXCENT,3,3),DERRAD(MXSP,MXCENT,3)
      COMMON /PCM_SURPOT/ POTCAVNUC(MXTS), POTCAVELE(MXTS)
