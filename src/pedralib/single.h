#if defined (VAR_SINGLE)
C
C    This file is for compiling on a NEC SX-4, and corresponds to CRAY
C    equivalences file.
C    It replaces all procedure calls to DOUBLE precision with the correspondig
C    SINGLE precision (but 64bit) routines.
C    This file will be included if the cpp flag (VAR_SINGLE) is defined, in all
C    files, except the files in pdpack (but also in linextra).
C
C
C    Let's hope that none of this names is used somewhere else!!!
C
C                                         Gilbert Hangartner, 24. July 1997

#define DASUM  SASUM
#define DAXPY  SAXPY
#define DCOPY  SCOPY
#define DDOT   SDOT
#define DNRM2  SNRM2
#define DSCAL  SSCAL
#define DSWAP  SSWAP
#define DGEMM  SGEMM
#define DGEMV  SGEMV
#define IDAMAX ISAMAX

C    Other equivalences, for routines found in the original pdpack/dsp.F and
C    pdpack/dge.F (inbetween moved to pdpack/linextra.F)
C    These routines are coming from linpack / eispack

#define DGECO  SGECO
#define DGEDI  SGEDI
#define DGEFA  SGEFA
#define DGESL  SGESL

#define DSPCO  SSPCO
#define DSPDI  SSPDI
#define DSPFA  SSPFA
#define DSPSL  SSPSL

#endif
