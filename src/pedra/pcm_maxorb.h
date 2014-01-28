!
!     File: maxorb.h
!
!     MXSHEL = maximum number of shells (insert shell definition here).
!     MXPRIM = maximum number of primitives.
!     MXCORB = maximum number of orbitals (possibly contracted).
!     MAXOCC = maximum number of occupied orbitals
!
!     IF you change any of these parameters you should do a "make depend"
!     and then rebuild the program using the command "make".
!
      INTEGER MXSHEL, MXPRIM, MXCORB, MXORBT, MAXOCC
      PARAMETER (MXSHEL = 750, MXPRIM = 8000, MXCORB = 1200,              &
                MAXOCC = 400, MXORBT = MXCORB*(MXCORB + 1)/2)
