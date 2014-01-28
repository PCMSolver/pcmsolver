!
!     file: mxcent.h
!
!     MXCENT = max number of nuclei + point charges + ghost orbital centers
!
!     IF you change MXCENT you should do a "make depend"
!     and then rebuild the program using the command "make".
!
!     In case of QM3 MXNEW is used to allocate memory in herrdn.F.
!     To run a QM3 calculation in most cases MXCENT will 
!     have to be around 2000 - 3000!!! Remember to set MXQM3 = MXCENT in
!     qm3.h!!!
!
      INTEGER MXNEW, MXCENT, MXCOOR
      PARAMETER (MXNEW =120, MXCENT = 120, MXCOOR = 3*MXCENT)
