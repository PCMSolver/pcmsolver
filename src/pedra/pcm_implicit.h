#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (VAR_REAL) || defined (SYS_T90)
      IMPLICIT REAL (A-H,O-Z)
#else
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#endif
