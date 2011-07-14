C ---- include/infpri.h ---
C      INFo on PRInt in SIRIUS module (primarily)
      LOGICAL P4FLAG,P6FLAG
      PARAMETER (NPFLAG = 60)
      COMMON /INFPRI/ IPRI4, IPRI6, IPRSIR,
     &                IPRCNO,IPRDIA,IPRSIG,IPRDNS,IPRSOL,
     &                IPRKAP,IPRMUL,IPRCIX,IPRRHF,IPRAVE,IPRFCK,
     &                MPRI4, MPRI6, MPRSIR,LIM_POPPRI,
     &                P4FLAG(NPFLAG),P6FLAG(NPFLAG)
C ---- end of include/infpri.h
