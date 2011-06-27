#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (SYS_T90)
CDIR$ IVDEP
#endif
#if defined (SYS_ALLIANT)
CVD$ NODEPCHK
#endif
#if defined (SYS_HAL)
C$DIR NO_RECURRENCE
#endif
#if defined (SYS_IBM)
Cvectorization note: ignore vector dependence
C*VDIR: IGNORE RECRDEPS
#endif
#if defined (SYS_NEC)
*VDIR NODEP
#endif
#if !defined (SYS_CRAY) && !defined (SYS_ALLIANT) && !defined (SYS_IBM) && !defined (SYS_T3D) && !defined (SYS_NEC) && !defined (SYS_HAL) && !defined (SYS_T90)
Cvectorization note: ignore vector dependence
#endif
