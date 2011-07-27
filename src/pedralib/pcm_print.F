C
C...   Dalton, Release 2.0, pdpack/printpkg.F
C...
C...   These routines are in the public domain and can be
C...   used freely in other programs.
C...
C
      SUBROUTINE OUTPAK_ (AMATRX,NROW,NCTL,LVPRI)
C.......................................................................
C Revised 16-Dec-1983 by Hans Jorgen Aa. Jensen.
C         16-Jun-1986 hjaaj ( removed Hollerith )
C
C OUTPAK_ PRINTS A REAL SYMMETRIC MATRIX STORED IN ROW-PACKED LOWER
C
C TRIANGULAR FORM (SEE DIAGRAM BELOW) IN FORMATTED FORM WITH NUMBERED
C
C ROWS AND COLUMNS.  THE INPUT IS AS FOLLOWS:
C
C        AMATRX(*)...........PACKED MATRIX
C
C        NROW................NUMBER OF ROWS TO BE OUTPUT_
C
C        NCTL................CARRIAGE CONTROL FLAG: 1 FOR SINGLE SPACE,
C                                                   2 FOR DOUBLE SPACE,
C                                                   3 FOR TRIPLE SPACE.
C
C THE MATRIX ELEMENTS ARE ARRANGED IN STORAGE AS FOLLOWS:
C
C        1
C        2    3
C        4    5    6
C        7    8    9   10
C       11   12   13   14   15
C       16   17   18   19   20   21
C       22   23   24   25   26   27   28
C       AND SO ON.
C
C OUTPAK_ IS SET UP TO HANDLE 6 COLUMNS/PAGE WITH A 6F20.14 FORMAT
C FOR THE COLUMNS.  IF A DIFFERENT NUMBER OF COLUMNS IS REQUIRED, CHANGE
C FORMATS 1000 AND 2000, AND INITIALIZE KCOL WITH THE NEW NUMBER OF
C COLUMNS.
C
C AUTHOR:  NELSON H.F. BEEBE, QUANTUM THEORY PROJECT, UNIVERSITY OF
C          FLORIDA, GAINESVILLE
C..........VERSION = 09/05/73/03
C.......................................................................
C
#include <pcm_implicit.h>
      DIMENSION AMATRX(*)
      INTEGER BEGIN
      CHARACTER*1 ASA(3),BLANK,CTL
      CHARACTER   PFMT*20, COLUMN*8
      PARAMETER (ZERO=0.D00, KCOLP=4, KCOLN=6)
      PARAMETER (FFMIN=1.D-3, FFMAX = 1.D3)
      DATA COLUMN/'Column  '/, ASA/' ', '0', '-'/, BLANK/' '/
C
      IF (NCTL .LT. 0) THEN
         KCOL = KCOLN
      ELSE
         KCOL = KCOLP
      END IF
      MCTL = ABS(NCTL)
      IF ((MCTL.LE.3).AND.(MCTL.GT.0)) THEN
         CTL = ASA(MCTL)
      ELSE
         CTL = BLANK
      END IF
C
      J = NROW*(NROW+1)/2
      AMAX = ZERO
      DO 5 I=1,J
         AMAX = MAX( AMAX, ABS(AMATRX(I)) )
    5 CONTINUE
      IF (AMAX .EQ. ZERO) THEN
         WRITE (LVPRI,'(/T6,A)') 'Zero matrix.'
         GO TO 200
      END IF
      IF (FFMIN .LE. AMAX .AND. AMAX .LE. FFMAX) THEN
C        use F output format
         PFMT = '(A1,I7,2X,8F15.8)'
      ELSE
C        use 1PD output format
         PFMT = '(A1,I7,2X,1P,8D15.6)'
      END IF
C
C LAST IS THE LAST COLUMN NUMBER IN THE ROW CURRENTLY BEING PRINTED
C
      LAST = MIN(NROW,KCOL)
C
C BEGIN IS THE FIRST COLUMN NUMBER IN THE ROW CURRENTLY BEING PRINTED.
C
C.....BEGIN NON STANDARD DO LOOP.
      BEGIN= 1
 1050 NCOL = 1
         WRITE (LVPRI,1000) (COLUMN,I,I = BEGIN,LAST)
         DO 40 K = BEGIN,NROW
            KTOTAL = (K*(K-1))/2 + BEGIN - 1
            DO 10 I = 1,NCOL
               IF (AMATRX(KTOTAL+I) .NE. ZERO) GO TO 20
   10       CONTINUE
            GO TO 30
   20       WRITE (LVPRI,PFMT) CTL,K,(AMATRX(J+KTOTAL),J=1,NCOL)
   30       IF (K .LT. (BEGIN+KCOL-1)) NCOL = NCOL + 1
   40    CONTINUE
         LAST = MIN(LAST+KCOL,NROW)
         BEGIN= BEGIN + NCOL
      IF (BEGIN.LE.NROW) GO TO 1050
  200 CONTINUE
      RETURN
C
 1000 FORMAT (/12X,6(3X,A6,I4,2X),(3X,A6,I4))
C2000 FORMAT (A1,'Row',I4,2X,1P,8D15.6)
C2000 FORMAT (A1,I7,2X,1P,8D15.6)
      END
C  /* Deck outpkb */
      SUBROUTINE OUTPKB_ (AMATRX,NDIM,NBLK,NCTL,LVPRI)
C.......................................................................
C
C OUTPKB_ is OUTPAK_ modified for blocking as specified by
C        NDIM(NBLK).  16-Dec-1983 Hans Jorgen Aa. Jensen.
C
C 16-Jun-1986 hjaaj ( removed Hollerith )
C
C OUTPAK_ PRINTS A REAL OMSYMMETRIC MATRIX STORED IN ROW-PACKED LOWER
C TRIANGULAR FORM (SEE DIAGRAM BELOW) IN FORMATTED FORM WITH NUMBERED
C ROWS AND COLUMNS.  THE INPUT IS AS FOLLOWS:
C
C        AMATRX(')...........PACKED MATRIX, blocked
C
C        NDIM(').............dimension of each block
C
C        NBLK................number of blocks
C
C        NCTL................CARRIAGE CONTROL FLAG: 1 FOR SINGLE SPACE,
C                                                   2 FOR DOUBLE SPACE,
C                                                   3 FOR TRIPLE SPACE.
C
C THE MATRIX ELEMENTS in a block ARE ARRANGED IN STORAGE AS
C FOLLOWS:
C
C        1
C        2    3
C        4    5    6
C        7    8    9   10
C       11   12   13   14   15
C       16   17   18   19   20   21
C       22   23   24   25   26   27   28
C       AND SO ON.
C
C OUTPAK_ IS SET UP TO HANDLE 6 COLUMNS/PAGE WITH A 6F20.14 FORMAT
C FOR THE COLUMNS.  IF A DIFFERENT NUMBER OF COLUMNS IS REQUIRED, CHANGE
C FORMATS 1000 AND 2000, AND INITIALIZE KCOL WITH THE NEW NUMBER OF
C COLUMNS.
C
C AUTHOR:  NELSON H.F. BEEBE, QUANTUM THEORY PROJECT, UNIVERSITY OF
C          FLORIDA, GAINESVILLE
C..........OUTPAK_ VERSION = 09/05/73/03
C.......................................................................
C
#include <pcm_implicit.h>
      INTEGER NDIM(NBLK),BEGIN
      DIMENSION AMATRX(*)
      CHARACTER*1 ASA(3),BLANK,CTL
      CHARACTER   PFMT*20, COLUMN*8
      PARAMETER (ZERO = 0.0D0, KCOLP=4, KCOLN=6)
      PARAMETER (FFMIN=1.D-3, FFMAX = 1.D3)
      DATA COLUMN/'Column  '/, ASA/' ', '0', '-'/, BLANK/' '/
C
      IF (NCTL .LT. 0) THEN
         KCOL = KCOLN
      ELSE
         KCOL = KCOLP
      END IF
      MCTL = ABS(NCTL)
      IF ((MCTL.LE.3).AND.(MCTL.GT.0)) THEN
         CTL = ASA(MCTL)
      ELSE
         CTL = BLANK
      END IF
C
      MATLN = 0
      DO 200 IBLK = 1,NBLK
         MATLN = MATLN + NDIM(IBLK)*(NDIM(IBLK)+1)/2
  200 CONTINUE
C
      AMAX = ZERO
      DO 205 I=1,MATLN
         AMAX = MAX( AMAX, ABS(AMATRX(I)) )
  205 CONTINUE
      IF (AMAX .EQ. ZERO) THEN
         WRITE (LVPRI,3000) NBLK
         GO TO 800
      END IF
      IF (FFMIN .LE. AMAX .AND. AMAX .LE. FFMAX) THEN
C        use F output format
         PFMT = '(A1,I7,2X,8F15.8)'
      ELSE
C        use 1PD output format
         PFMT = '(A1,I7,2X,1P,8D15.6)'
      END IF
C
      IOFF = 0
      DO 100 IBLK = 1,NBLK
         IDIM = NDIM(IBLK)
      IF (IDIM.EQ.0) GO TO 100
         IIDIM = IDIM*(IDIM+1)/2
C
         DO 5 I=1,IIDIM
            IF (AMATRX(IOFF+I).NE.ZERO) GO TO 15
    5    CONTINUE
         WRITE (LVPRI,1100) IBLK
         GO TO 90
   15    CONTINUE
         WRITE (LVPRI,1200) IBLK
C
C LAST IS THE LAST COLUMN NUMBER IN THE ROW CURRENTLY BEING PRINTED
C
         LAST = MIN(IDIM,KCOL)
C
C BEGIN IS THE FIRST COLUMN NUMBER IN THE ROW CURRENTLY BEING PRINTED.
C
         BEGIN = 1
C.....BEGIN NON STANDARD DO LOOP.
 1050       NCOL = 1
            WRITE (LVPRI,1000) (COLUMN,I,I = BEGIN,LAST)
            KTOTAL = IOFF + BEGIN*(BEGIN+1)/2 - 1
            DO 40 K = BEGIN,IDIM
               DO 10 I = 1,NCOL
                  IF (AMATRX(KTOTAL+I) .NE. ZERO) GO TO 20
   10          CONTINUE
               GO TO 30
   20          WRITE (LVPRI,PFMT) CTL,K,(AMATRX(KTOTAL+J),J=1,NCOL)
   30          IF (K .LT. (BEGIN+KCOL-1)) NCOL = NCOL + 1
               KTOTAL = KTOTAL + K
   40       CONTINUE
            LAST = MIN(LAST+KCOL,IDIM)
            BEGIN=BEGIN+NCOL
         IF (BEGIN.LE.IDIM) GO TO 1050
C
   90    IOFF = IOFF + IIDIM
  100 CONTINUE
C
  800 CONTINUE
      RETURN
 3000 FORMAT (/5X,'All',I3,' blocks zero matrices.')
 1100 FORMAT (/5X,'*** Block',I3,' zero matrix. ***')
 1200 FORMAT (/5X,'*** Block',I3,' ***')
 1000 FORMAT (/12X,6(3X,A6,I4,2X),(3X,A6,I4))
C2000 FORMAT (A1,'Row',I4,2X,1P,8D15.6)
C2000 FORMAT (A1,I7,2X,1P,8D15.6)
      END
C  /* Deck outpks */
      SUBROUTINE OUTPKS_ (AMATRX,NDIM,NBLK,IREPO,NCTL,LVPRI)
C.......................................................................
C
C OUTPKS_ is OUTPAK_ modified for blocking as specified by
C        NDIM(NBLK).  16-Dec-1983 Hans Jorgen Aa. Jensen.
C        14-Dec-1988 tuh , accepts not totally symmetric elements
C
C 16-Jun-1986 hjaaj ( removed Hollerith )
C
C OUTPAK_ PRINTS A REAL OMSYMMETRIC MATRIX STORED IN ROW-PACKED LOWER
C TRIANGULAR FORM (SEE DIAGRAM BELOW) IN FORMATTED FORM WITH NUMBERED
C ROWS AND COLUMNS.  THE INPUT IS AS FOLLOWS:
C
C        AMATRX(')...........PACKED MATRIX, blocked
C
C        NDIM(').............dimension of each block
C
C        NBLK................number of blocks
C
C        NCTL................CARRIAGE CONTROL FLAG: 1 FOR SINGLE SPACE,
C                                                   2 FOR DOUBLE SPACE,
C                                                   3 FOR TRIPLE SPACE.
C
C THE MATRIX ELEMENTS in a block ARE ARRANGED IN STORAGE AS
C FOLLOWS:
C
C        1
C        2    3
C        4    5    6
C        7    8    9   10
C       11   12   13   14   15
C       16   17   18   19   20   21
C       22   23   24   25   26   27   28
C       AND SO ON.
C
C OUTPAK_ IS SET UP TO HANDLE 6 COLUMNS/PAGE WITH A 6F20.14 FORMAT
C FOR THE COLUMNS.  IF A DIFFERENT NUMBER OF COLUMNS IS REQUIRED, CHANGE
C FORMATS 1000 AND 2000, AND INITIALIZE KCOL WITH THE NEW NUMBER OF
C COLUMNS.
C
C AUTHOR:  NELSON H.F. BEEBE, QUANTUM THEORY PROJECT, UNIVERSITY OF
C          FLORIDA, GAINESVILLE
C..........OUTPAK_ VERSION = 09/05/73/03
C.......................................................................
C
#include <pcm_implicit.h>
      INTEGER NDIM(NBLK),BEGIN
      DIMENSION AMATRX(*)
      CHARACTER*1 ASA(3),BLANK,CTL
      CHARACTER   PFMT*20, COLUMN*8
      PARAMETER (D0 = 0.0D0, KCOLP=4, KCOLN=6)
      PARAMETER (FFMIN=1.D-3, FFMAX = 1.D3)
      DATA COLUMN/'Column  '/, ASA/' ', '0', '-'/, BLANK/' '/
#include <pcm_ibtfun.h>
C
      IF (IREPO .NE. 1) THEN
         IBLK = 1
         DO 500 IREPB = 0, NBLK - 1
            IREPA = IBTXOR(IREPO - 1,IREPB)
            IF (IREPA .GT. IREPB) THEN
               NDIMA = NDIM(IREPA + 1)
               NDIMB = NDIM(IREPB + 1)
               WRITE (LVPRI,'(/A,2I5)') ' Symmetries: ',IREPA+1,IREPB+1
               WRITE (LVPRI,'(/A,2I5)') ' Dimensions: ',NDIMA,NDIMB
               IF (NDIMA.GT.0 .AND. NDIMB.GT.0) THEN
                  CALL OUTPUT_(AMATRX(IBLK),1,NDIMA,1,NDIMB,
     *                        NDIMA,NDIMB,NCTL,LVPRI)
                  IBLK = IBLK + NDIMA*NDIMB
               END IF
            ENDIF
 500     CONTINUE
         GO TO 800
      END IF
      IF (NCTL .LT. 0) THEN
         KCOL = KCOLN
      ELSE
         KCOL = KCOLP
      END IF
      MCTL = ABS(NCTL)
      IF ((MCTL.LE.3).AND.(MCTL.GT.0)) THEN
         CTL = ASA(MCTL)
      ELSE
         CTL = BLANK
      END IF
C
      MATLN = 0
      DO 200 IBLK = 1,NBLK
         MATLN = MATLN + NDIM(IBLK)*(NDIM(IBLK)+1)/2
  200 CONTINUE
C
      AMAX = D0
      DO 205 I=1,MATLN
         AMAX = MAX( AMAX, ABS(AMATRX(I)) )
  205 CONTINUE
      IF (AMAX .EQ. D0) THEN
         WRITE (LVPRI,3000) NBLK
         GO TO 800
      END IF
      IF (FFMIN .LE. AMAX .AND. AMAX .LE. FFMAX) THEN
C        use F output format
         PFMT = '(A1,I7,2X,8F15.8)'
      ELSE
C        use 1PD output format
         PFMT = '(A1,I7,2X,1P,8D15.6)'
      END IF
C
      IOFF = 0
      DO 100 IBLK = 1,NBLK
         IDIM = NDIM(IBLK)
      IF (IDIM.EQ.0) GO TO 100
         IIDIM = IDIM*(IDIM+1)/2
C
         DO 5 I=1,IIDIM
            IF (AMATRX(IOFF+I).NE.D0) GO TO 15
    5    CONTINUE
         WRITE (LVPRI,1100) IBLK
         GO TO 90
   15    CONTINUE
         WRITE (LVPRI,1200) IBLK
C
C LAST IS THE LAST COLUMN NUMBER IN THE ROW CURRENTLY BEING PRINTED
C
         LAST = MIN(IDIM,KCOL)
C
C BEGIN IS THE FIRST COLUMN NUMBER IN THE ROW CURRENTLY BEING PRINTED.
C
         BEGIN = 1
C.....BEGIN NON STANDARD DO LOOP.
 1050       NCOL = 1
            WRITE (LVPRI,1000) (COLUMN,I,I = BEGIN,LAST)
            KTOTAL = IOFF + BEGIN*(BEGIN+1)/2 - 1
            DO 40 K = BEGIN,IDIM
               DO 10 I = 1,NCOL
                  IF (AMATRX(KTOTAL+I) .NE. D0) GO TO 20
   10          CONTINUE
               GO TO 30
   20          WRITE (LVPRI,PFMT) CTL,K,(AMATRX(KTOTAL+J),J=1,NCOL)
   30          IF (K .LT. (BEGIN+KCOL-1)) NCOL = NCOL + 1
               KTOTAL = KTOTAL + K
   40       CONTINUE
            LAST = MIN(LAST+KCOL,IDIM)
            BEGIN=BEGIN+NCOL
         IF (BEGIN.LE.IDIM) GO TO 1050
C
   90    IOFF = IOFF + IIDIM
  100 CONTINUE
C
  800 CONTINUE
      RETURN
 3000 FORMAT (/5X,'All',I3,' blocks zero matrices.')
 1100 FORMAT (/5X,'*** Block',I3,' zero matrix. ***')
 1200 FORMAT (/5X,'*** Block',I3,' ***')
 1000 FORMAT (/12X,6(3X,A6,I4,2X),(3X,A6,I4))
      END
C  /* Deck outptb */
      SUBROUTINE OUTPTB_ (AMATRX,NDIM,NBLK,NCTL,LVPRI)
C.......................................................................
C
C OUTPTB_ is OUTPUT_ modified for blocking as specified by
C        NDIM(NBLK).  16-Dec-1983 Hans Jorgen Aa. Jensen.
C
C 16-Jun-1986 hjaaj ( removed Hollerith )
C
C OUTPUT_ PRINTS A REAL MATRIX IN FORMATTED FORM WITH NUMBERED ROWS
C AND COLUMNS.  THE INPUT IS AS FOLLOWS;
C
C        AMATRX(',').........MATRIX TO BE OUTPUT_, blocked
C
C        NDIM(').............dimension of each block
C
C        NBLK................number of blocks
C
C        NCTL................CARRIAGE CONTROL FLAG; 1 FOR SINGLE SPACE
C                                                   2 FOR DOUBLE SPACE
C                                                   3 FOR TRIPLE SPACE
C
C THE PARAMETERS THAT FOLLOW MATRIX ARE ALL OF TYPE INTEGER*4.  THE
C PROGRAM IS SET UP TO HANDLE 5 COLUMNS/PAGE WITH A 1P,5D24.15 FORMAT
C FOR 1THE COLUMNS.  IF A DIFFERENT NUMBER OF COLUMNS IS REQUIRED,
C CHANGE FORMATS 1000 AND 2000, AND INITIALIZE KCOL WITH THE NEW NUMBER
C OF COLUMNS.
C
C AUTHOR;  NELSON H.F. BEEBE, QUANTUM THEORY PROJECT, UNIVERSITY OF
C          FLORIDA, GAINESVILLE
C REVISED; FEBRUARY 26, 1971
C
C.......................................................................
C
#include <pcm_implicit.h>
      DIMENSION AMATRX(*)
      INTEGER NDIM(NBLK), BEGIN, KCOL
      CHARACTER*1 ASA(3), BLANK, CTL
      CHARACTER   PFMT*20, COLUMN*8
      PARAMETER (ZERO=0.D00, KCOLP=4, KCOLN=6)
      PARAMETER (FFMIN=1.D-3, FFMAX = 1.D3)
      DATA COLUMN/'Column  '/, ASA/' ', '0', '-'/, BLANK/' '/
C
      IF (NCTL .LT. 0) THEN
         KCOL = KCOLN
      ELSE
         KCOL = KCOLP
      END IF
      MCTL = ABS(NCTL)
      IF ((MCTL.LE.3).AND.(MCTL.GT.0)) THEN
         CTL = ASA(MCTL)
      ELSE
         CTL = BLANK
      END IF
C
      MATDIM = 0
      DO 200 IBLK = 1,NBLK
         MATDIM = MATDIM + NDIM(IBLK)*NDIM(IBLK)
  200 CONTINUE
C
      AMAX = ZERO
      DO 205 I=1,MATDIM
         AMAX = MAX( AMAX, ABS(AMATRX(I)) )
  205 CONTINUE
      IF (AMAX .EQ. ZERO) THEN
         WRITE (LVPRI,3000) NBLK
         GO TO 800
      END IF
      IF (FFMIN .LE. AMAX .AND. AMAX .LE. FFMAX) THEN
C        use F output format
         PFMT = '(A1,I7,2X,8F15.8)'
      ELSE
C        use 1PD output format
         PFMT = '(A1,I7,2X,1P,8D15.6)'
      END IF
C
      IOFF = 0
        DO 100 IBLK = 1,NBLK
        IDIM = NDIM(IBLK)
        IF (IDIM.EQ.0) GO TO 100
        I2DIM = IDIM*IDIM
          DO 10 I=1,I2DIM
          IF (AMATRX(IOFF+I).NE.ZERO) GO TO 15
   10     CONTINUE
        WRITE (LVPRI,1100) IBLK
        GO TO 90
   15   CONTINUE
        WRITE (LVPRI,1200) IBLK
        LAST = MIN(IDIM,KCOL)
          DO 2 BEGIN = 1,IDIM,KCOL
          WRITE (LVPRI,1000) (COLUMN,I,I = BEGIN,LAST)
            DO 1 K = 1,IDIM
            IOFFK = IOFF + K
              DO 4 I=BEGIN,LAST
              IF (AMATRX(IOFFK+(I-1)*IDIM).NE.ZERO) GO TO 5
    4         CONTINUE
            GO TO 1
    5       WRITE (LVPRI,PFMT) CTL,K,(AMATRX(IOFFK+(I-1)*IDIM),
     *        I = BEGIN,LAST)
    1       CONTINUE
    2     LAST = MIN(LAST+KCOL,IDIM)
   90   IOFF = IOFF + I2DIM
  100   CONTINUE
C
  800 CONTINUE
      RETURN
 3000 FORMAT (/5X,'All',I3,' blocks zero matrices.')
 1100 FORMAT (/5X,'*** Block',I3,' zero matrix. ***')
 1200 FORMAT (/5X,'*** Block',I3,' ***')
 1000 FORMAT (/12X,6(3X,A6,I4,2X),(3X,A6,I4))
C2000 FORMAT (A1,'Row',I4,2X,1P,8D15.6)
C2000 FORMAT (A1,I7,2X,1P,8D15.6)
      END
C  /* Deck output */
      SUBROUTINE OUTPUT_(AMATRX,ROWLOW,ROWHI,COLLOW,COLHI,ROWDIM,COLDIM,
     *                   NCTL,LVPRI)
C.......................................................................
C Revised 15-Dec-1983 by Hans Jorgen Aa. Jensen.
C         16-Jun-1986 hjaaj ( removed Hollerith )
C
C OUTPUT_ PRINTS A REAL MATRIX IN FORMATTED FORM WITH NUMBERED ROWS
C AND COLUMNS.  THE INPUT IS AS FOLLOWS;
C
C        AMATRX(',').........MATRIX TO BE OUTPUT_
C
C        ROWLOW..............ROW NUMBER AT WHICH OUTPUT_ IS TO BEGIN
C
C        ROWHI...............ROW NUMBER AT WHICH OUTPUT_ IS TO END
C
C        COLLOW..............COLUMN NUMBER AT WHICH OUTPUT_ IS TO BEGIN
C
C        COLHI...............COLUMN NUMBER AT WHICH OUTPUT_ IS TO END
C
C        ROWDIM..............ROW DIMENSION OF AMATRX(',')
C
C        COLDIM..............COLUMN DIMENSION OF AMATRX(',')
C
C        NCTL................CARRIAGE CONTROL FLAG; 1 FOR SINGLE SPACE
C                                                   2 FOR DOUBLE SPACE
C                                                   3 FOR TRIPLE SPACE
C                            hjaaj: negative for 132 col width
C
C THE PARAMETERS THAT FOLLOW MATRIX ARE ALL OF TYPE INTEGER*4.  THE
C PROGRAM IS SET UP TO HANDLE 5 COLUMNS/PAGE WITH A 1P,5D24.15 FORMAT
C FOR THE COLUMNS.  IF A DIFFERENT NUMBER OF COLUMNS IS REQUIRED,
C CHANGE FORMATS 1000 AND 2000, AND INITIALIZE KCOL WITH THE NEW NUMBER
C OF COLUMNS.
C
C AUTHOR;  NELSON H.F. BEEBE, QUANTUM THEORY PROJECT, UNIVERSITY OF
C          FLORIDA, GAINESVILLE
C REVISED; FEBRUARY 26, 1971
C
C.......................................................................
C
#include <pcm_implicit.h>
      INTEGER   ROWLOW,ROWHI,COLLOW,COLHI,ROWDIM,COLDIM,BEGIN,KCOL
      DIMENSION AMATRX(ROWDIM,COLDIM)
      CHARACTER*1 ASA(3), BLANK, CTL
      CHARACTER   PFMT*20, COLUMN*8
      PARAMETER (ZERO=0.D00, KCOLP=5, KCOLN=8)
      PARAMETER (FFMIN=1.D-3, FFMAX = 1.D3)
      DATA COLUMN/'Column  '/, BLANK/' '/, ASA/' ', '0', '-'/
C
      IF (ROWHI.LT.ROWLOW) GO TO 3
      IF (COLHI.LT.COLLOW) GO TO 3
C
      AMAX = ZERO
      DO 10 J = COLLOW,COLHI
         DO 10 I = ROWLOW,ROWHI
            AMAX = MAX( AMAX, ABS(AMATRX(I,J)) )
   10 CONTINUE
      IF (AMAX .EQ. ZERO) THEN
         WRITE (LVPRI,'(/T6,A)') 'Zero matrix.'
         GO TO 3
      END IF
      IF (FFMIN .LE. AMAX .AND. AMAX .LT. FFMAX) THEN
C        use F output format
         PFMT = '(A1,I7,2X,8F14.8)'
         thrpri = 0.5D-8
      ELSE
C        use 1PD output format
         PFMT = '(A1,I7,2X,1P,8D14.6)'
         thrpri = 1.0D-8*AMAX
      END IF
C
      IF (NCTL .LT. 0) THEN
         KCOL = KCOLN
      ELSE
         KCOL = KCOLP
      END IF
      MCTL = ABS(NCTL)
      IF ((MCTL.LE.3).AND.(MCTL.GT.0)) THEN
         CTL = ASA(MCTL)
      ELSE
         CTL = BLANK
      END IF
C
      LAST = MIN(COLHI,COLLOW+KCOL-1)
      DO 2 BEGIN = COLLOW,COLHI,KCOL
         WRITE (LVPRI,1000) (COLUMN,I,I = BEGIN,LAST)
         DO 1 K = ROWLOW,ROWHI
            DO 4 I = BEGIN,LAST
               IF (abs(AMATRX(K,I)).gt.thrpri) GO TO 5
    4       CONTINUE
         GO TO 1
    5       WRITE (LVPRI,PFMT) CTL,K,(AMATRX(K,I), I = BEGIN,LAST)
    1    CONTINUE
    2 LAST = MIN(LAST+KCOL,COLHI)
    3 RETURN
 1000 FORMAT (/10X,8(4X,A6,I4))
C2000 FORMAT (A1,'Row',I4,2X,1P,8D14.6)
C2000 FORMAT (A1,I7,2X,1P,8D14.6)
      END
      SUBROUTINE OUTPUT_x(AMATRX,ROWLOW,ROWHI,COLLOW,COLHI,ROWDIM,
     *                   COLDIM,NCTL,LVPRI)
C.......................................................................
C Revised 15-Dec-1983 by Hans Jorgen Aa. Jensen.
C         16-Jun-1986 hjaaj ( removed Hollerith )
C OUTPUT_x based on OUTPUT 2007, Hans Jorgen Aa. Jensen
C         prints 7F10.3/D10.3 instead of 5F14.6/D14.6
C
C OUTPUT_ PRINTS A REAL MATRIX IN FORMATTED FORM WITH NUMBERED ROWS
C AND COLUMNS.  THE INPUT IS AS FOLLOWS;
C
C        AMATRX(',').........MATRIX TO BE OUTPUT_
C
C        ROWLOW..............ROW NUMBER AT WHICH OUTPUT_ IS TO BEGIN
C
C        ROWHI...............ROW NUMBER AT WHICH OUTPUT_ IS TO END
C
C        COLLOW..............COLUMN NUMBER AT WHICH OUTPUT_ IS TO BEGIN
C
C        COLHI...............COLUMN NUMBER AT WHICH OUTPUT_ IS TO END
C
C        ROWDIM..............ROW DIMENSION OF AMATRX(',')
C
C        COLDIM..............COLUMN DIMENSION OF AMATRX(',')
C
C        NCTL................CARRIAGE CONTROL FLAG; 1 FOR SINGLE SPACE
C                                                   2 FOR DOUBLE SPACE
C                                                   3 FOR TRIPLE SPACE
C                            hjaaj: negative for 132 col width
C
C THE PARAMETERS THAT FOLLOW MATRIX ARE ALL OF TYPE INTEGER*4.  THE
C PROGRAM IS SET UP TO HANDLE 5 COLUMNS/PAGE WITH A 1P,5D24.15 FORMAT
C FOR THE COLUMNS.  IF A DIFFERENT NUMBER OF COLUMNS IS REQUIRED,
C CHANGE FORMATS 1000 AND 2000, AND INITIALIZE KCOL WITH THE NEW NUMBER
C OF COLUMNS.
C
C AUTHOR;  NELSON H.F. BEEBE, QUANTUM THEORY PROJECT, UNIVERSITY OF
C          FLORIDA, GAINESVILLE
C REVISED; FEBRUARY 26, 1971
C
C.......................................................................
C
#include <pcm_implicit.h>
      INTEGER   ROWLOW,ROWHI,COLLOW,COLHI,ROWDIM,COLDIM,BEGIN,KCOL
      DIMENSION AMATRX(ROWDIM,COLDIM)
      CHARACTER*1 ASA(3), BLANK, CTL
      CHARACTER   PFMT*20, COLUMN*8
      PARAMETER (ZERO=0.D00, KCOLP=7, KCOLN=10)
      PARAMETER (FFMIN=1.D-3, FFMAX = 1.D3)
      DATA COLUMN/'Column  '/, BLANK/' '/, ASA/' ', '0', '-'/
C
      IF (ROWHI.LT.ROWLOW) GO TO 3
      IF (COLHI.LT.COLLOW) GO TO 3
C
      AMAX = ZERO
      DO 10 J = COLLOW,COLHI
         DO 10 I = ROWLOW,ROWHI
            AMAX = MAX( AMAX, ABS(AMATRX(I,J)) )
   10 CONTINUE
      IF (AMAX .EQ. ZERO) THEN
         WRITE (LVPRI,'(/T6,A)') 'Zero matrix.'
         GO TO 3
      END IF
      IF (FFMIN .LE. AMAX .AND. AMAX .LE. FFMAX) THEN
C        use F output format
         PFMT = '(A1,I5,2X,10F10.3)'
      ELSE
C        use 1PD output format
         PFMT = '(A1,I5,2X,1P,10D10.3)'
      END IF
C
      IF (NCTL .LT. 0) THEN
         KCOL = KCOLN
      ELSE
         KCOL = KCOLP
      END IF
      MCTL = ABS(NCTL)
      IF ((MCTL.LE.3).AND.(MCTL.GT.0)) THEN
         CTL = ASA(MCTL)
      ELSE
         CTL = BLANK
      END IF
C
      LAST = MIN(COLHI,COLLOW+KCOL-1)
      DO 2 BEGIN = COLLOW,COLHI,KCOL
         WRITE (LVPRI,1000) (COLUMN,I,I = BEGIN,LAST)
         DO 1 K = ROWLOW,ROWHI
            DO 4 I = BEGIN,LAST
               IF (abs(AMATRX(K,I)).GT.0.5D-3) GO TO 5
    4       CONTINUE
         GO TO 1
    5       WRITE (LVPRI,PFMT) CTL,K,(AMATRX(K,I), I = BEGIN,LAST)
    1    CONTINUE
    2 LAST = MIN(LAST+KCOL,COLHI)
    3 RETURN
 1000 FORMAT (/9X,10(A6,I3,1X),(3X,A6,I4))
C2000 FORMAT (A1,'Row',I4,2X,1P,8D15.6)
C2000 FORMAT (A1,I7,2X,1P,8D15.6)
      END
C  /* Deck prtab */
      SUBROUTINE PRTAB_(NTABLE,TABLE,TEXT,LVPRI)
C
C 28-Dec-1989 Hans Joergen Aa. Jensen
C
C Purpose: print tables of text (e.g. from input parsing routines)
C
      CHARACTER*(*) TABLE(NTABLE), TEXT
      LTEXT = LEN(TEXT)
      LTEXT = MIN(70,LTEXT)
      WRITE (LVPRI,'(//1X,A/1X,70A1/)') TEXT(1:LTEXT),('=',I=1,LTEXT)
      DO 100 I = 1,NTABLE
         IF (INDEX(TABLE(I),'XXXX') .EQ. 0) THEN
            WRITE (LVPRI,'(T6,A)') TABLE(I)
         END IF
  100 CONTINUE
      WRITE (LVPRI,'()')
      RETURN
      END
C  /* Deck prvec */
      SUBROUTINE PRVEC_(NDIM,VEC,INCVEC,THRESH,MAXLIN,LUOUT)
C
C 19-Aug-1989 Hans Joergen Aa. Jensen
C
C   print VEC(1:NDIM*INCVEC:INCVEC)
C
C   NDIM      : Number of elements in vector VEC
C   VEC(:)    : Vector to be printed
C   INCVEC    : Increment between each element in vector VEC
C   THRESH    : Print threshold for vector with unit norm
C               (if THRESH .lt. 0 then -THRESH is used without
C                renormalization).
C   MAXLIN    : max. lines of output with vector elements
C   LUOUT     : output unit
C
C
#include <pcm_implicit.h>
      DIMENSION VEC(*)
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0)
      PARAMETER (D1LOW = D1 - 1.0D-10, D1HGH = D1 + 1.0D-10)
      LOGICAL   VSCALE
C
C     Test input
C
      NERR = 0
      IF (NDIM   .LE. 0) THEN
         WRITE (LUOUT,'(/5X,A)')
     &     'No print from PRVEC_ because NDIM .le. 0'
         NERR = NERR + 1
      END IF
      IF (INCVEC .LE. 0) THEN
         WRITE (LUOUT,'(/5X,A)')
     &      'No print from PRVEC_ because INCVEC .le. 0'
         NERR = NERR + 1
      END IF
      IF (NERR .GT. 0) RETURN
C
C
C
      C2NRM = DNORM2__(NDIM,VEC,INCVEC)
      IF (THRESH .LE. D0) THEN
         THRES2 = -THRESH
         VSCALE = .FALSE.
      ELSE
         THRES2 =  THRESH * C2NRM
         VSCALE = (C2NRM .LT. D1LOW .OR. C2NRM .GT. D1HGH)
      END IF
      IF (VSCALE) THEN
         WRITE (LUOUT,'(/2A,1P,D12.4)')
     &      ' Print of vector elements (vector scaled to unit norm)',
     &      ' larger than',ABS(THRESH)
         SCALE = D1 / C2NRM
      ELSE
         WRITE (LUOUT,'(/A,1P,D12.4)')
     &      ' Print of vector elements larger than',ABS(THRESH)
         SCALE = D1
      END IF
C
      IPR   = 0
      IZER  = 0
      C2SUM = D0
C
      DO 300 I = 1, NDIM
         NA = 1 + (I-1)*INCVEC
         IF (ABS(VEC(NA)).LE.THRES2 .OR. IPR .GE. MAXLIN) THEN
            C2SUM = C2SUM + VEC(NA)*VEC(NA)
            IF (VEC(NA) .EQ. D0) IZER = IZER + 1
         ELSE
            IF (MOD(IPR,5) .EQ. 0) WRITE (LUOUT,'()')
            IPR = IPR + 1
            IF (VSCALE) THEN
               WRITE(LUOUT,50)I,VEC(NA),SCALE*VEC(NA)
            ELSE
               WRITE(LUOUT,60)I,VEC(NA)
            END IF
   50       FORMAT(3X,'Element',I10,3X,'coefficient',1P,D10.2,0P,
     &             3X,'scaled to unit norm',F10.6)
   60       FORMAT(6X,'Element',I12,3X,'coefficient',F20.10)
         END IF
  300 CONTINUE
      IF (IPR .GE. MAXLIN) THEN
C
C     *** We have reached the print limit
C
         WRITE (LUOUT,910) IPR
      END IF
  910 FORMAT(/' Print limit of',I6,' elements has been reached.')
      C2SUM = SQRT(C2SUM)
      IF (IZER .EQ. NDIM) THEN
         WRITE (LUOUT,920) NDIM
      ELSE
         WRITE (LUOUT,930) NDIM,NDIM-IPR,IZER,C2SUM,C2NRM
      END IF
  920 FORMAT(/' Length of vector                      :',I10,
     *       /' All elements are zero.',/)
  930 FORMAT(/' Length of vector                      :',I10,
     *       /' Number of elements not printed        :',I10,
     *       /' Number of zero elements               :',I10,
     *       /' Total norm of coefficients not printed:',F10.6,
     *       /' (the coefficients are normalized to    ',F10.6,')',/)
      RETURN
C
C End of PRVEC_
C
      END
C  /* Deck prmgn */
      SUBROUTINE PRMGN_(NDIM,VEC,INCVEC,NPOT,LUOUT)
C
C 19-Aug-1989 hjaaj
C
C   NDIM      : Number of elements in vector VEC
C   VEC(:)    : Vector to be printed
C   INCVEC    : Increment between each element in vector VEC
C   NPOT      : IBASE**-NPOT is lowest magnitude considered
C   LUOUT     : output unit
C
C
#include <pcm_implicit.h>
      DIMENSION VEC(*)
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0)
      PARAMETER (IBASE = 10)
C     IBASE     : Base for magnitude
C
      FACMGN = IBASE
      FACMGN = D1 / FACMGN
C
C     Test input
C
      NERR = 0
      IF (NDIM   .LE. 0) THEN
         WRITE (LUOUT,*) 'No print from PRMGN_ because NDIM .le. 0'
         NERR = NERR + 1
      END IF
      IF (INCVEC .LE. 0) THEN
         WRITE (LUOUT,*) 'No print from PRMGN_ because INCVEC .le. 0'
         NERR = NERR + 1
      END IF
      IF (NERR .GT. 0) RETURN
C
C
C
      C2NRM = DNORM2__(NDIM,VEC,INCVEC)
      IZER  = 0
      NLAST = NDIM*INCVEC
      DO 100 NA = 1, NLAST, INCVEC
         IF (VEC(NA) .EQ. D0) IZER = IZER + 1
  100 CONTINUE
C
C.... SIZE OF COEFFICIENTS
C
      WRITE(LUOUT,'(/A)')  '  Magnitude of coefficients '
      WRITE(LUOUT,'( A)')  ' ==========================='
      IF (IZER .EQ. NDIM) THEN
         WRITE(LUOUT,'(/A)') ' All elements are zero.'
         GO TO 9999
      END IF
C
      IF (C2NRM .GT. (D1 + 1.0 D-10) ) THEN
         WRITE (LUOUT,'(/A,1P,E23.12,A,/A)')
     &      ' (Vector norm is',C2NRM,')',
     &      ' (Range is scaled to a vector norm of 1)'
      END IF
      XMIN  = MAX(C2NRM,D1)
C     ... max possible coefficient is C2NRM
C
      WRITE (LUOUT,'(/2A,/2A)')
     &   '     Range        # of elements',
     &   '     Norm squared         Accumulated',
     &   ' -------------   --------------',
     &   '   ---------------     ---------------'
C
C Output will look like:
C
C     Range        # of elements     Norm squared         Accumulated
C -------------   --------------   ---------------     ---------------
C 10- 1 to  1                2      1.12345678E-01      1.12345678E-01
C
      C2NRM = D0
      IDET  = IZER
C. LOOP OVER INTER_VALS
      IPOT  = 0
  200 CONTINUE
         CLNORM = D0
         INRANG = 0
         XMAX   = XMIN
         IF (IPOT .EQ. 0) XMAX = 2*XMAX
C        ... to make certain we don't omit any because of round-off
         XMIN   = XMIN * FACMGN
         DO 300 NA = 1, NLAST, INCVEC
           IF( ABS(VEC(NA)) .LE. XMAX  .AND.
     &         ABS(VEC(NA)) .GT. XMIN ) THEN
                 INRANG = INRANG + 1
                 CLNORM = CLNORM + VEC(NA) ** 2
           END IF
  300    CONTINUE
         C2NRM = C2NRM + CLNORM
C
C
         IF (INRANG .GT. 0) THEN
            IF (IPOT .EQ. 0) THEN
               WRITE(LUOUT,'(I3,A,I14,1P,2E20.8)')
     &            IBASE,'- 1 to  1   ',INRANG,CLNORM,C2NRM
            ELSE
               WRITE(LUOUT,'(I3,A,I2,A,I3,A,I2,I14,1P,2E20.8)')
     &            IBASE,'-',IPOT+1,' to',IBASE,'-',IPOT,
     &            INRANG,CLNORM,C2NRM
            END IF
         END IF
C
C
         IDET = IDET + INRANG
         IPOT = IPOT + 1
      IF( IDET .LT. NDIM .AND. IPOT .LT. NPOT ) GOTO 200
C
      ISML = NDIM - IDET
      IF (ISML .GT. 0) WRITE (LUOUT,'(A,I13)') ' Other non-zero ',ISML
      IF (IZER .GT. 0) WRITE (LUOUT,'(A,I13)') ' Exact zero     ',IZER
                     WRITE (LUOUT,'(/A,I13/)') ' Total          ',NDIM
C
C
 9999 CONTINUE
      RETURN
C
C End of PRMGN_
C
      END
C  /* Deck pnzvec */
      FUNCTION PNZVEC_(NDIM,VEC,INCVEC,THRZER)
C
C Copyright 24-Mar-1993 Hans Joergen Aagaard Jensen
C   return % non-zero elements of VEC(1:NDIM*INCVEC:INCVEC)
C
C   NDIM      : Number of elements in vector VEC
C   VEC(:)    : Vector
C   INCVEC    : Increment between each element in vector VEC
C   THRZER    : Threshold for zero elements
C
C
#include <pcm_implicit.h>
      DIMENSION VEC(*)
#include <pcm_priunit.h>
      PARAMETER (D0 = 0.0D0, D100 = 100.0D0)
C
C     Test input
C
      NERR = 0
      IF (NDIM   .LE. 0) THEN
         WRITE (LVPRI,'(/5X,A,I10)')
     &     'No percentage from PNZVEC_ because NDIM .le. 0; NDIM =',NDIM
         NERR = NERR + 1
      END IF
      IF (INCVEC .LE. 0) THEN
         WRITE (LVPRI,'(/5X,A,I10)')
     &      'No percentage from PNZVEC_ because INCVEC .le. 0;INCVEC =',
     &      INCVEC
         NERR = NERR + 1
      END IF
      IF (NERR .GT. 0) THEN
         PNZVEC_ = -D100
         RETURN
      END IF
C
      THRESH = MAX(D0,THRZER)
      NNZ = 0
      DO 100 I = 1,NDIM*INCVEC,INCVEC
         IF (ABS(VEC(I)) .GT. THRESH) NNZ = NNZ + 1
  100 CONTINUE
      DNZ    = NNZ
      DNDIM  = NDIM
      PNZVEC_ = DNZ*D100 / DNDIM
      RETURN
      END
C  /* Deck coutput */
      SUBROUTINE COUTPUT__(AMATRX,ROWLOW,ROWHI,COLLOW,COLHI,ROWDIM,
     *                   COLDIM,NCTL,LVPRI)
C.......................................................................
C Revised 15-Dec-1983 by Hans Jorgen Aa. Jensen.
C         16-Jun-1986 hjaaj ( removed Hollerith )
C         20-Mar-2001 ov    complex matrix verison
C
C OUTPUT_ PRINTS A REAL MATRIX IN FORMATTED FORM WITH NUMBERED ROWS
C AND COLUMNS.  THE INPUT IS AS FOLLOWS;
C
C        AMATRX(',').........MATRIX TO BE OUTPUT_
C
C        ROWLOW..............ROW NUMBER AT WHICH OUTPUT_ IS TO BEGIN
C
C        ROWHI...............ROW NUMBER AT WHICH OUTPUT_ IS TO END
C
C        COLLOW..............COLUMN NUMBER AT WHICH OUTPUT_ IS TO BEGIN
C
C        COLHI...............COLUMN NUMBER AT WHICH OUTPUT_ IS TO END
C
C        ROWDIM..............ROW DIMENSION OF AMATRX(',')
C
C        COLDIM..............COLUMN DIMENSION OF AMATRX(',')
C
C        NCTL................CARRIAGE CONTROL FLAG; 1 FOR SINGLE SPACE
C                                                   2 FOR DOUBLE SPACE
C                                                   3 FOR TRIPLE SPACE
C
C THE PARAMETERS THAT FOLLOW MATRIX ARE ALL OF TYPE INTEGER*4.  THE
C PROGRAM IS SET UP TO HANDLE 5 COLUMNS/PAGE WITH A 1P,5D24.15 FORMAT
C FOR THE COLUMNS.  IF A DIFFERENT NUMBER OF COLUMNS IS REQUIRED,
C CHANGE FORMATS 1000 AND 2000, AND INITIALIZE KCOL WITH THE NEW NUMBER
C OF COLUMNS.
C
C AUTHOR;  NELSON H.F. BEEBE, QUANTUM THEORY PROJECT, UNIVERSITY OF
C          FLORIDA, GAINESVILLE
C REVISED; FEBRUARY 26, 1971
C
C.......................................................................
C
      IMPLICIT DOUBLE COMPLEX (A-H,O-Z)
      INTEGER   ROWLOW,ROWHI,COLLOW,COLHI,ROWDIM,COLDIM,BEGIN,KCOL
      DIMENSION AMATRX(ROWDIM,COLDIM)
      CHARACTER*1 ASA(3), BLANK, CTL
      CHARACTER   PFMT*20, COLUMN*8
      DOUBLE PRECISION ZERO, FFMIN, FFMAX
      PARAMETER (ZERO=0.D00, KCOLP=4, KCOLN=6)
      PARAMETER (FFMIN=1.D-3, FFMAX = 1.D3)
      DATA COLUMN/'Column  '/, BLANK/' '/, ASA/' ', '0', '-'/
      DOUBLE PRECISION AMAX
C
      IF (ROWHI.LT.ROWLOW) GO TO 3
      IF (COLHI.LT.COLLOW) GO TO 3
C
      AMAX = ZERO
      DO 10 J = COLLOW,COLHI
         DO 10 I = ROWLOW,ROWHI
            AMAX = MAX( AMAX, ABS(AMATRX(I,J)) )
   10 CONTINUE
      IF (AMAX .EQ. ZERO) THEN
         WRITE (LVPRI,'(/T6,A)') 'Zero matrix.'
         GO TO 3
      END IF
      IF (FFMIN .LE. AMAX .AND. AMAX .LE. FFMAX) THEN
C        use F output format
         PFMT = '(A1,I3,2X,4(2F10.6))'
      ELSE
C        use 1PD output format
         PFMT = '(A1,I7,2X,1P,4(2D12.6))'
      END IF
C
      IF (NCTL .LT. 0) THEN
         KCOL = KCOLN
      ELSE
         KCOL = KCOLP
      END IF
      MCTL = ABS(NCTL)
      IF ((MCTL.LE.3).AND.(MCTL.GT.0)) THEN
         CTL = ASA(MCTL)
      ELSE
         CTL = BLANK
      END IF
C
      LAST = MIN(COLHI,COLLOW+KCOL-1)
      DO 2 BEGIN = COLLOW,COLHI,KCOL
         WRITE (LVPRI,1000) (COLUMN,I,I = BEGIN,LAST)
         DO 1 K = ROWLOW,ROWHI
            DO 4 I = BEGIN,LAST
               IF (AMATRX(K,I).NE.ZERO) GO TO 5
    4       CONTINUE
         GO TO 1
    5       WRITE (LVPRI,PFMT) CTL,K,(AMATRX(K,I), I = BEGIN,LAST)
    1    CONTINUE
    2 LAST = MIN(LAST+KCOL,COLHI)
    3 RETURN
 1000 FORMAT (/5X,6(3X,A6,I4,7X),(3X,A6,I4))
C2000 FORMAT (A1,'Row',I4,2X,1P,8D15.6)
C2000 FORMAT (A1,I7,2X,1P,8D15.6)
      END
C --- end of printpkg.F ---
