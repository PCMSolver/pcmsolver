C
C     File: infind.h - orbital index vectors for sirius module + abacus and response;
C                      vectors generally defined in sirius/sirset.F
C
C     IROW      : for indexing of triangular packed matrices
C     ISMO/ISAO : symmetry of each MO/AO
C     ISW       : reorders to inactive-active-secondary ordering from symmetry ordering
C     ISX       : reorders back
C     ICH       : if neg, inactive no.; if pos., active no.; if zero, empty orbital
C     IOBTYP    : define type of each MO with JTFRO, JTINAC, JTACT, or JTSEC
C     NSM       : symmetry of active orbitals
C     IACTYP    : RAS block for each active orbital
C     ISSMO     : super symmetry of each MO (defined in sirius/sirave.F)
C     ISSORD    : reorder array to super symmetry ordering (defined in sirius/sirave.F)
C
      INTEGER LIROW, JTFRO, JTINAC, JTACT, JTSEC, JTFRFR, JTINFR,
     &        JTACFR, JTSEFR, JTININ, JTACIN, JTACAC, JTSEIN, JTSEAC,
     &        JTSESE
      PARAMETER (LIROW  = 3*MXCORB + 3)
C     Definition of readable codes for orbital types and 2-index orbital types
      PARAMETER (JTFRO  = 1, JTINAC = 2, JTACT  = 3, JTSEC  = 4)
      PARAMETER (JTFRFR =-1, JTINFR =-2, JTACFR =-3, JTSEFR =-4,
     &           JTININ = 1, JTACIN = 2, JTACAC = 3, JTSEIN = 4,
     &           JTSEAC = 5, JTSESE = 6)
C
      INTEGER IROW, ISMO, ISAO, ISW, ISX, ICH, IOBTYP,
     &        NSM, IACTYP, ISSMO, ISSORD
      COMMON /INFIND/ IROW(LIROW),  ISMO(MXCORB),ISAO(MXCORB),
     &                ISW(MXCORB),  ISX(MXCORB),
     &                ICH(MXCORB),  IOBTYP(MXCORB),
     &                NSM(MAXASH),  IACTYP(MAXASH),
     &                ISSMO(MXCORB),ISSORD(MXCORB)
