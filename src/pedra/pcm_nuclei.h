      INTEGER NUCPRE, NUCNUM, NUCDEG, ISTBNU, NCTOT,                &
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
      INTEGER MULBSI
      COMMON /PCM_MULBAS/ MULBSI(MXCENT)
!     MULBAS has been added for multiple basis sets (WK/UniKA/31-10-2002).
