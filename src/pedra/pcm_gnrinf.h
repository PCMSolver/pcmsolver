! File : gnrinf.h
!
!     -*- mode: fortran; fortran-continuation-string: "&" -*-
!     File: gnrinf.h -- general information for DALTON
!
      LOGICAL TESTIN, OPTWLK, RNHERM, RNSIRI, RNABAC, GEOCNV,       &
             HRINPC, SRINPC, RDINPC, RDMLIN, PARCAL, DIRCAL,        &
             WRINDX, WLKREJ, WALKIN, RNRESP, USRIPR, SEGBAS,        &
             DOCCSD, OPTNEW, NEWSYM, NEWBAS, NEWPRP, RELCAL,        &
             TOTSYM, NMWALK, DKTRAN, GEOALL, WESTA,  SEGAUX,        &
             ERFEXP, DOSRIN, SRINTS, CHI1ST, QM3, QMMM
!
      COMMON /PCM_GNRINF/ 
             ! double:                                              &
             GRADML, PANAS,  CHIVAL, THR_REDFAC,                    &
             ! integer:                                             &
             KCHARG, ITERNR, ITERMX, IPRUSR, LENBAS,                &
             ! logical:                                             &
             TESTIN, OPTWLK, RNHERM, RNSIRI, RNABAC, GEOCNV,        &
             HRINPC, SRINPC, RDINPC, RDMLIN, PARCAL, DIRCAL,        &
             WRINDX, WLKREJ, WALKIN, RNRESP, USRIPR, SEGBAS,        &
             DOCCSD, OPTNEW, NEWSYM, NEWBAS, NEWPRP, RELCAL,        &
             TOTSYM, NMWALK, DKTRAN, GEOALL, WESTA,  SEGAUX,        &
             ERFEXP, DOSRIN, SRINTS, CHI1ST, QM3, QMMM              &

      INTEGER LBASDIR
      PARAMETER (LBASDIR = 600)
      CHARACTER*(LBASDIR) BASDIR
      COMMON /PCM_GNRCHR/ BASDIR
! --- end of gnrinf.h ---
