#if defined(__CVERSION__)
struct common_priunit {
  int lucmd, lupri, luerr, luw4, ninfo, nwarn, iprerr, lupot;
};
extern struct common_priunit priunit_;
#else
C     FILE: priunit.h
      CHARACTER*80 SEPARATOR
      PARAMETER (SEPARATOR = '----------------------------------------'
     &                     //'----------------------------------------')
      INTEGER LUCMD, LUPRI, LUERR , LUW4, NINFO, NWARN, IPRERR, LUPOT
      COMMON /PRIUNIT/
     &        LUCMD, LUPRI, LUERR , LUW4, NINFO, NWARN, IPRERR, LUPOT
#endif
