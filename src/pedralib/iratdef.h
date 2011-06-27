C     FILE: iratdef.h
C     
C     IRAT  = (real word length) / (integer word length)
C     LRAT  = (real word length) / (logical word length)
C
      INTEGER IRAT, LRAT
#if defined (VAR_INT64)
C     using INTEGER*8 (64 bit integers) as default ...
      PARAMETER (IRAT = 1, LRAT = 1)
#else
      PARAMETER (IRAT = 2, LRAT = 2)
#endif
