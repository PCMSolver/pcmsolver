C
C     File: maxorb.h
C
C     MXSHEL = maximum number of shells (insert shell definition here).
C     MXPRIM = maximum number of primitives.
C     MXCORB = maximum number of orbitals (possibly contracted).
C     MAXOCC = maximum number of occupied orbitals
C
C     IF you change any of these parameters you should do a "make depend"
C     and then rebuild the program using the command "make".
C
      INTEGER MXSHEL, MXPRIM, MXCORB, MXORBT, MAXOCC
      PARAMETER (MXSHEL = 750, MXPRIM = 8000, MXCORB = 1200,
     *           MAXOCC = 400, MXORBT = MXCORB*(MXCORB + 1)/2)
