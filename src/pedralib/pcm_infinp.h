C File : infinp.h
C
c  -*- mode:fortran; fortran-continuation-string: "&" -*-
C     this common block contains general Sirius input read in sirius/sirinp.F
C     (specified under **WAVE FUNCTION).  /hjaaj Oct 2003
C     - SCF specific input is in scbrhf.h
C     - orbital specifications are in inforb.h
C
      INTEGER         NFLAG, MAXRTS, MXFELT
      PARAMETER (NFLAG = 80, MAXRTS = 100, MXFELT = 20)
C
      INTEGER         NFIELD, ISPIN,NMCAVER,ISTATE,LSYM,NACTEL, MCTYPE,
     *                LSOLMX,NLMSOL,NELMN1,NELMX1,NELMN3,NELMX3,
     *                LROOTS,NROOTS,IROOT,
     *                NOROT        ,IMOORD,
     *                IORTO,ICI0,KDEL,ICHECK,NTIT,
     *                MAXMAC,MAXMIC,MAXJT,MAXCIT,MAXUIT,MAXAPM,MAXABS,
     *                ITRLVL,ITRFIN,JCHSYM,JCHORB,
     *                NROOCI,ISTACI, MXCIMA, ICICNO,IMCCNO
      COMMON /PCM_INTINP/ NFIELD, ISPIN,ISTATE,LSYM,NACTEL, MCTYPE,
     *                LSOLMX,NLMSOL,NELMN1,NELMX1,NELMN3,NELMX3,
     *                LROOTS,NROOTS,IROOT(MAXRTS),
     *                NOROT(MXCORB),IMOORD(MXCORB),
     *                IORTO,ICI0,KDEL,ICHECK,NTIT,
     *                MAXMAC,MAXMIC,MAXJT,MAXCIT,MAXUIT,MAXAPM,MAXABS,
     *                ITRLVL,ITRFIN,JCHSYM,JCHORB,
     *                NROOCI,ISTACI, MXCIMA, ICICNO,IMCCNO, NMCAVER
C
      LOGICAL         FLAG,  DOSCF, DOMP2, DOCINO,DOCI,  DOMC,  DORSP,
     &                FCVORB,LNOROT,LMOORD,DIRFCK,CORHOL,CORRLX,RESPHP,
     &                JOLSEN,ABAIPH,INERSI,INERSF,DODFT, DONEVPT,HSROHF,
     &                BOYORB,PIPORB,ADDMP2
C     variables for srDFT /hjaaj
      LOGICAL         DOCISRDFT,DOHFSRDFT,DOMCSRDFT,ADDSRI,SRHYBR
      COMMON /PCM_LOGINP/FLAG(NFLAG),DOSCF,DOMP2,DOCINO,DOCI,DOMC,DORSP,
     &                FCVORB,LNOROT,LMOORD,DIRFCK,CORHOL,CORRLX,RESPHP,
     &                JOLSEN,ABAIPH,INERSI,INERSF,DODFT,DONEVPT,HSROHF,
     &                BOYORB,PIPORB,ADDMP2,
     &
     &                DOCISRDFT,DOHFSRDFT,DOMCSRDFT,ADDSRI,SRHYBR
      LOGICAL         SUPSYM, DORHF
      EQUIVALENCE (SUPSYM,FLAG(17)), (DOSCF,DORHF)
C
      REAL*8          SPIN, POTNUC, EPSOL,EPSTAT,EPPN,RSOL,
     &                THRGRD, THRPWF, THRCI, THRMC, THRCGR,
     &                EFIELD, TITMOL, CMAXMO, THROVL,
     &                THRSSY, DEFLVL, WEIGHT_MCAVER
C     variables for srDFT /hjaaj
      REAL*8          THRCIDFT
      COMMON /PCM_RELINP/ SPIN, POTNUC, EPSOL,EPSTAT,EPPN,RSOL(3),
     &                THRGRD, THRPWF, THRCI, THRMC, THRCGR,
     &                EFIELD(MXFELT), TITMOL(12,2), CMAXMO, THROVL,
     &                THRSSY, DEFLVL, WEIGHT_MCAVER(MAXRTS),
     &
     &                THRCIDFT
C
      CHARACTER*60 TITLE
      CHARACTER*4  CENT,   TYPE
      CHARACTER*8  LFIELD
      COMMON /PCM_CHRINP/ TITLE(6), CENT(MXCORB), TYPE(MXCORB),
     &                LFIELD(MXFELT)
C-- end of infinp.h --
