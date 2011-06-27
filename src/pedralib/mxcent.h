C
C     file: mxcent.h
C
C     MXCENT = max number of nuclei + point charges + ghost orbital centers
C
C     IF you change MXCENT you should do a "make depend"
C     and then rebuild the program using the command "make".
C
C     In case of QM3 MXNEW is used to allocate memory in herrdn.F.
C     To run a QM3 calculation in most cases MXCENT will 
C     have to be around 2000 - 3000!!! Remember to set MXQM3 = MXCENT in
C     qm3.h!!!
C
      INTEGER MXNEW, MXCENT, MXCOOR
      PARAMETER (MXNEW =120, MXCENT = 120, MXCOOR = 3*MXCENT)
