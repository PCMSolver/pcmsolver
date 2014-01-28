!     FILE: iratdef.h
!     
!     IRAT  = (real word length) / (integer word length)
!     LRAT  = (real word length) / (logical word length)
!
      INTEGER IRAT, LRAT
#if defined (VAR_INT64)
!     using INTEGER*8 (64 bit integers) as default ...
      PARAMETER (IRAT = 1, LRAT = 1)
#else
      PARAMETER (IRAT = 2, LRAT = 2)
#endif
