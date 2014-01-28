
!...   Dalton, Release 2.0, pdpack/printpkg.F
!...
!...   These routines are in the public domain and can be
!...   used freely in other programs.
!...

    SUBROUTINE OUTPAK_ (AMATRX,NROW,NCTL,LVPRI)
!.......................................................................
! Revised 16-Dec-1983 by Hans Jorgen Aa. Jensen.
!         16-Jun-1986 hjaaj ( removed Hollerith )

! OUTPAK_ PRINTS A REAL SYMMETRIC MATRIX STORED IN ROW-PACKED LOWER

! TRIANGULAR FORM (SEE DIAGRAM BELOW) IN FORMATTED FORM WITH NUMBERED

! ROWS AND COLUMNS.  THE INPUT IS AS FOLLOWS:

!        AMATRX(*)...........PACKED MATRIX

!        NROW................NUMBER OF ROWS TO BE OUTPUT_

!        NCTL................CARRIAGE CONTROL FLAG: 1 FOR SINGLE SPACE,
!                                                   2 FOR DOUBLE SPACE,
!                                                   3 FOR TRIPLE SPACE.

! THE MATRIX ELEMENTS ARE ARRANGED IN STORAGE AS FOLLOWS:

!        1
!        2    3
!        4    5    6
!        7    8    9   10
!       11   12   13   14   15
!       16   17   18   19   20   21
!       22   23   24   25   26   27   28
!       AND SO ON.

! OUTPAK_ IS SET UP TO HANDLE 6 COLUMNS/PAGE WITH A 6F20.14 FORMAT
! FOR THE COLUMNS.  IF A DIFFERENT NUMBER OF COLUMNS IS REQUIRED, CHANGE
! FORMATS 1000 AND 2000, AND INITIALIZE KCOL WITH THE NEW NUMBER OF
! COLUMNS.

! AUTHOR:  NELSON H.F. BEEBE, QUANTUM THEORY PROJECT, UNIVERSITY OF
!          FLORIDA, GAINESVILLE
!..........VERSION = 09/05/73/03
!.......................................................................

#include <pcm_implicit.h>
    DIMENSION AMATRX(*)
    INTEGER :: BEGIN
    CHARACTER(1) :: ASA(3),BLANK,CTL
    CHARACTER   PFMT*20, COLUMN*8
    PARAMETER (ZERO=0.D00, KCOLP=4, KCOLN=6)
    PARAMETER (FFMIN=1.D-3, FFMAX = 1.D3)
    DATA COLUMN/'Column  '/, ASA/' ', '0', '-'/, BLANK/' '/

    IF (NCTL < 0) THEN
        KCOL = KCOLN
    ELSE
        KCOL = KCOLP
    END IF
    MCTL = ABS(NCTL)
    IF ((MCTL <= 3) .AND. (MCTL > 0)) THEN
        CTL = ASA(MCTL)
    ELSE
        CTL = BLANK
    END IF

    J = NROW*(NROW+1)/2
    AMAX = ZERO
    DO 5 I=1,J
        AMAX = MAX( AMAX, ABS(AMATRX(I)) )
    5 END DO
    IF (AMAX == ZERO) THEN
        WRITE (LVPRI,'(/T6,A)') 'Zero matrix.'
        GO TO 200
    END IF
    IF (FFMIN <= AMAX .AND. AMAX <= FFMAX) THEN
    !        use F output format
        PFMT = '(A1,I7,2X,8F15.8)'
    ELSE
    !        use 1PD output format
        PFMT = '(A1,I7,2X,1P,8D15.6)'
    END IF

! LAST IS THE LAST COLUMN NUMBER IN THE ROW CURRENTLY BEING PRINTED

    LAST = MIN(NROW,KCOL)

! BEGIN IS THE FIRST COLUMN NUMBER IN THE ROW CURRENTLY BEING PRINTED.

!.....BEGIN NON STANDARD DO LOOP.
    BEGIN= 1
    1050 NCOL = 1
    WRITE (LVPRI,1000) (COLUMN,I,I = BEGIN,LAST)
    DO 40 K = BEGIN,NROW
        KTOTAL = (K*(K-1))/2 + BEGIN - 1
        DO 10 I = 1,NCOL
            IF (AMATRX(KTOTAL+I) /= ZERO) GO TO 20
        10 END DO
        GO TO 30
    20 WRITE (LVPRI,PFMT) CTL,K,(AMATRX(J+KTOTAL),J=1,NCOL)
    30 IF (K < (BEGIN+KCOL-1)) NCOL = NCOL + 1
    40 END DO
    LAST = MIN(LAST+KCOL,NROW)
    BEGIN= BEGIN + NCOL
    IF (BEGIN <= NROW) GO TO 1050
    200 CONTINUE
    RETURN

    1000 FORMAT (/12X,6(3X,A6,I4,2X),(3X,A6,I4))
!2000 FORMAT (A1,'Row',I4,2X,1P,8D15.6)
!2000 FORMAT (A1,I7,2X,1P,8D15.6)
    END SUBROUTINE OUTPAK_
!  /* Deck outpkb */
    SUBROUTINE OUTPKB_ (AMATRX,NDIM,NBLK,NCTL,LVPRI)
!.......................................................................

! OUTPKB_ is OUTPAK_ modified for blocking as specified by
!        NDIM(NBLK).  16-Dec-1983 Hans Jorgen Aa. Jensen.

! 16-Jun-1986 hjaaj ( removed Hollerith )

! OUTPAK_ PRINTS A REAL OMSYMMETRIC MATRIX STORED IN ROW-PACKED LOWER
! TRIANGULAR FORM (SEE DIAGRAM BELOW) IN FORMATTED FORM WITH NUMBERED
! ROWS AND COLUMNS.  THE INPUT IS AS FOLLOWS:

!        AMATRX(')...........PACKED MATRIX, blocked

!        NDIM(').............dimension of each block

!        NBLK................number of blocks

!        NCTL................CARRIAGE CONTROL FLAG: 1 FOR SINGLE SPACE,
!                                                   2 FOR DOUBLE SPACE,
!                                                   3 FOR TRIPLE SPACE.

! THE MATRIX ELEMENTS in a block ARE ARRANGED IN STORAGE AS
! FOLLOWS:

!        1
!        2    3
!        4    5    6
!        7    8    9   10
!       11   12   13   14   15
!       16   17   18   19   20   21
!       22   23   24   25   26   27   28
!       AND SO ON.

! OUTPAK_ IS SET UP TO HANDLE 6 COLUMNS/PAGE WITH A 6F20.14 FORMAT
! FOR THE COLUMNS.  IF A DIFFERENT NUMBER OF COLUMNS IS REQUIRED, CHANGE
! FORMATS 1000 AND 2000, AND INITIALIZE KCOL WITH THE NEW NUMBER OF
! COLUMNS.

! AUTHOR:  NELSON H.F. BEEBE, QUANTUM THEORY PROJECT, UNIVERSITY OF
!          FLORIDA, GAINESVILLE
!..........OUTPAK_ VERSION = 09/05/73/03
!.......................................................................

#include <pcm_implicit.h>
    INTEGER :: NDIM(NBLK),BEGIN
    DIMENSION AMATRX(*)
    CHARACTER(1) :: ASA(3),BLANK,CTL
    CHARACTER   PFMT*20, COLUMN*8
    PARAMETER (ZERO = 0.0D0, KCOLP=4, KCOLN=6)
    PARAMETER (FFMIN=1.D-3, FFMAX = 1.D3)
    DATA COLUMN/'Column  '/, ASA/' ', '0', '-'/, BLANK/' '/

    IF (NCTL < 0) THEN
        KCOL = KCOLN
    ELSE
        KCOL = KCOLP
    END IF
    MCTL = ABS(NCTL)
    IF ((MCTL <= 3) .AND. (MCTL > 0)) THEN
        CTL = ASA(MCTL)
    ELSE
        CTL = BLANK
    END IF

    MATLN = 0
    DO 200 IBLK = 1,NBLK
        MATLN = MATLN + NDIM(IBLK)*(NDIM(IBLK)+1)/2
    200 END DO

    AMAX = ZERO
    DO 205 I=1,MATLN
        AMAX = MAX( AMAX, ABS(AMATRX(I)) )
    205 END DO
    IF (AMAX == ZERO) THEN
        WRITE (LVPRI,3000) NBLK
        GO TO 800
    END IF
    IF (FFMIN <= AMAX .AND. AMAX <= FFMAX) THEN
    !        use F output format
        PFMT = '(A1,I7,2X,8F15.8)'
    ELSE
    !        use 1PD output format
        PFMT = '(A1,I7,2X,1P,8D15.6)'
    END IF

    IOFF = 0
    DO 100 IBLK = 1,NBLK
        IDIM = NDIM(IBLK)
        IF (IDIM == 0) GO TO 100
        IIDIM = IDIM*(IDIM+1)/2
            
        DO 5 I=1,IIDIM
            IF (AMATRX(IOFF+I) /= ZERO) GO TO 15
        5 END DO
    WRITE (LVPRI,1100) IBLK
    GO TO 90
    15 CONTINUE
    WRITE (LVPRI,1200) IBLK
        
! LAST IS THE LAST COLUMN NUMBER IN THE ROW CURRENTLY BEING PRINTED
        
    LAST = MIN(IDIM,KCOL)
        
! BEGIN IS THE FIRST COLUMN NUMBER IN THE ROW CURRENTLY BEING PRINTED.
        
    BEGIN = 1
!.....BEGIN NON STANDARD DO LOOP.
    1050 NCOL = 1
    WRITE (LVPRI,1000) (COLUMN,I,I = BEGIN,LAST)
    KTOTAL = IOFF + BEGIN*(BEGIN+1)/2 - 1
    DO 40 K = BEGIN,IDIM
        DO 10 I = 1,NCOL
            IF (AMATRX(KTOTAL+I) /= ZERO) GO TO 20
        10 END DO
    GO TO 30
    20 WRITE (LVPRI,PFMT) CTL,K,(AMATRX(KTOTAL+J),J=1,NCOL)
    30 IF (K < (BEGIN+KCOL-1)) NCOL = NCOL + 1
    KTOTAL = KTOTAL + K
    40 END DO
    LAST = MIN(LAST+KCOL,IDIM)
    BEGIN=BEGIN+NCOL
    IF (BEGIN <= IDIM) GO TO 1050
        
    90 IOFF = IOFF + IIDIM
    100 END DO

    800 CONTINUE
    RETURN
    3000 FORMAT (/5X,'All',I3,' blocks zero matrices.')
    1100 FORMAT (/5X,'*** Block',I3,' zero matrix. ***')
    1200 FORMAT (/5X,'*** Block',I3,' ***')
    1000 FORMAT (/12X,6(3X,A6,I4,2X),(3X,A6,I4))
! 2000 FORMAT (A1,'Row',I4,2X,1P,8D15.6)
! 2000 FORMAT (A1,I7,2X,1P,8D15.6)
    END SUBROUTINE OUTPKB_
!  /* Deck outpks */
    SUBROUTINE OUTPKS_ (AMATRX,NDIM,NBLK,IREPO,NCTL,LVPRI)
!.......................................................................

! OUTPKS_ is OUTPAK_ modified for blocking as specified by
!        NDIM(NBLK).  16-Dec-1983 Hans Jorgen Aa. Jensen.
!        14-Dec-1988 tuh , accepts not totally symmetric elements

! 16-Jun-1986 hjaaj ( removed Hollerith )

! OUTPAK_ PRINTS A REAL OMSYMMETRIC MATRIX STORED IN ROW-PACKED LOWER
! TRIANGULAR FORM (SEE DIAGRAM BELOW) IN FORMATTED FORM WITH NUMBERED
! ROWS AND COLUMNS.  THE INPUT IS AS FOLLOWS:

!        AMATRX(')...........PACKED MATRIX, blocked

!        NDIM(').............dimension of each block

!        NBLK................number of blocks

!        NCTL................CARRIAGE CONTROL FLAG: 1 FOR SINGLE SPACE,
!                                                   2 FOR DOUBLE SPACE,
!                                                   3 FOR TRIPLE SPACE.

! THE MATRIX ELEMENTS in a block ARE ARRANGED IN STORAGE AS
! FOLLOWS:

!        1
!        2    3
!        4    5    6
!        7    8    9   10
!       11   12   13   14   15
!       16   17   18   19   20   21
!       22   23   24   25   26   27   28
!       AND SO ON.

! OUTPAK_ IS SET UP TO HANDLE 6 COLUMNS/PAGE WITH A 6F20.14 FORMAT
! FOR THE COLUMNS.  IF A DIFFERENT NUMBER OF COLUMNS IS REQUIRED, CHANGE
! FORMATS 1000 AND 2000, AND INITIALIZE KCOL WITH THE NEW NUMBER OF
! COLUMNS.

! AUTHOR:  NELSON H.F. BEEBE, QUANTUM THEORY PROJECT, UNIVERSITY OF
!          FLORIDA, GAINESVILLE
!..........OUTPAK_ VERSION = 09/05/73/03
!.......................................................................

#include <pcm_implicit.h>
    INTEGER :: NDIM(NBLK),BEGIN
    DIMENSION AMATRX(*)
    CHARACTER(1) :: ASA(3),BLANK,CTL
    CHARACTER   PFMT*20, COLUMN*8
    PARAMETER (D0 = 0.0D0, KCOLP=4, KCOLN=6)
    PARAMETER (FFMIN=1.D-3, FFMAX = 1.D3)
    DATA COLUMN/'Column  '/, ASA/' ', '0', '-'/, BLANK/' '/
#include <pcm_ibtfun.h>

    IF (IREPO /= 1) THEN
        IBLK = 1
        DO 500 IREPB = 0, NBLK - 1
            IREPA = IBTXOR(IREPO - 1,IREPB)
            IF (IREPA > IREPB) THEN
                NDIMA = NDIM(IREPA + 1)
                NDIMB = NDIM(IREPB + 1)
                WRITE (LVPRI,'(/A,2I5)') ' Symmetries: ',IREPA+1,IREPB+1
                WRITE (LVPRI,'(/A,2I5)') ' Dimensions: ',NDIMA,NDIMB
                IF (NDIMA > 0 .AND. NDIMB > 0) THEN
                    CALL OUTPUT_(AMATRX(IBLK),1,NDIMA,1,NDIMB, &
                    NDIMA,NDIMB,NCTL,LVPRI)
                    IBLK = IBLK + NDIMA*NDIMB
                END IF
            ENDIF
        500 END DO
    GO TO 800
    END IF
    IF (NCTL < 0) THEN
        KCOL = KCOLN
    ELSE
        KCOL = KCOLP
    END IF
    MCTL = ABS(NCTL)
    IF ((MCTL <= 3) .AND. (MCTL > 0)) THEN
        CTL = ASA(MCTL)
    ELSE
        CTL = BLANK
    END IF

    MATLN = 0
    DO 200 IBLK = 1,NBLK
        MATLN = MATLN + NDIM(IBLK)*(NDIM(IBLK)+1)/2
    200 END DO

    AMAX = D0
    DO 205 I=1,MATLN
        AMAX = MAX( AMAX, ABS(AMATRX(I)) )
    205 END DO
    IF (AMAX == D0) THEN
        WRITE (LVPRI,3000) NBLK
        GO TO 800
    END IF
    IF (FFMIN <= AMAX .AND. AMAX <= FFMAX) THEN
    !        use F output format
        PFMT = '(A1,I7,2X,8F15.8)'
    ELSE
    !        use 1PD output format
        PFMT = '(A1,I7,2X,1P,8D15.6)'
    END IF

    IOFF = 0
    DO 100 IBLK = 1,NBLK
        IDIM = NDIM(IBLK)
        IF (IDIM == 0) GO TO 100
        IIDIM = IDIM*(IDIM+1)/2
            
        DO 5 I=1,IIDIM
            IF (AMATRX(IOFF+I) /= D0) GO TO 15
        5 END DO
    WRITE (LVPRI,1100) IBLK
    GO TO 90
    15 CONTINUE
    WRITE (LVPRI,1200) IBLK
        
! LAST IS THE LAST COLUMN NUMBER IN THE ROW CURRENTLY BEING PRINTED
        
    LAST = MIN(IDIM,KCOL)
        
! BEGIN IS THE FIRST COLUMN NUMBER IN THE ROW CURRENTLY BEING PRINTED.
        
    BEGIN = 1
!.....BEGIN NON STANDARD DO LOOP.
    1050 NCOL = 1
    WRITE (LVPRI,1000) (COLUMN,I,I = BEGIN,LAST)
    KTOTAL = IOFF + BEGIN*(BEGIN+1)/2 - 1
    DO 40 K = BEGIN,IDIM
        DO 10 I = 1,NCOL
            IF (AMATRX(KTOTAL+I) /= D0) GO TO 20
        10 END DO
    GO TO 30
    20 WRITE (LVPRI,PFMT) CTL,K,(AMATRX(KTOTAL+J),J=1,NCOL)
    30 IF (K < (BEGIN+KCOL-1)) NCOL = NCOL + 1
    KTOTAL = KTOTAL + K
    40 END DO
    LAST = MIN(LAST+KCOL,IDIM)
    BEGIN=BEGIN+NCOL
    IF (BEGIN <= IDIM) GO TO 1050
        
    90 IOFF = IOFF + IIDIM
    100 END DO

    800 CONTINUE
    RETURN
    3000 FORMAT (/5X,'All',I3,' blocks zero matrices.')
    1100 FORMAT (/5X,'*** Block',I3,' zero matrix. ***')
    1200 FORMAT (/5X,'*** Block',I3,' ***')
    1000 FORMAT (/12X,6(3X,A6,I4,2X),(3X,A6,I4))
    END SUBROUTINE OUTPKS_
!  /* Deck outptb */
    SUBROUTINE OUTPTB_ (AMATRX,NDIM,NBLK,NCTL,LVPRI)
!.......................................................................

! OUTPTB_ is OUTPUT_ modified for blocking as specified by
!        NDIM(NBLK).  16-Dec-1983 Hans Jorgen Aa. Jensen.

! 16-Jun-1986 hjaaj ( removed Hollerith )

! OUTPUT_ PRINTS A REAL MATRIX IN FORMATTED FORM WITH NUMBERED ROWS
! AND COLUMNS.  THE INPUT IS AS FOLLOWS;

!        AMATRX(',').........MATRIX TO BE OUTPUT_, blocked

!        NDIM(').............dimension of each block

!        NBLK................number of blocks

!        NCTL................CARRIAGE CONTROL FLAG; 1 FOR SINGLE SPACE
!                                                   2 FOR DOUBLE SPACE
!                                                   3 FOR TRIPLE SPACE

! THE PARAMETERS THAT FOLLOW MATRIX ARE ALL OF TYPE INTEGER*4.  THE
! PROGRAM IS SET UP TO HANDLE 5 COLUMNS/PAGE WITH A 1P,5D24.15 FORMAT
! FOR 1THE COLUMNS.  IF A DIFFERENT NUMBER OF COLUMNS IS REQUIRED,
! CHANGE FORMATS 1000 AND 2000, AND INITIALIZE KCOL WITH THE NEW NUMBER
! OF COLUMNS.

! AUTHOR;  NELSON H.F. BEEBE, QUANTUM THEORY PROJECT, UNIVERSITY OF
!          FLORIDA, GAINESVILLE
! REVISED; FEBRUARY 26, 1971

!.......................................................................

#include <pcm_implicit.h>
    DIMENSION AMATRX(*)
    INTEGER :: NDIM(NBLK), BEGIN, KCOL
    CHARACTER(1) :: ASA(3), BLANK, CTL
    CHARACTER   PFMT*20, COLUMN*8
    PARAMETER (ZERO=0.0D00, KCOLP=4, KCOLN=6)
    PARAMETER (FFMIN=1.0D-03, FFMAX = 1.0D03)
    DATA COLUMN/'Column  '/, ASA/' ', '0', '-'/, BLANK/' '/

    IF (NCTL < 0) THEN
        KCOL = KCOLN
    ELSE
        KCOL = KCOLP
    END IF
    MCTL = ABS(NCTL)
    IF ((MCTL <= 3) .AND. (MCTL > 0)) THEN
        CTL = ASA(MCTL)
    ELSE
        CTL = BLANK
    END IF

    MATDIM = 0
    DO 200 IBLK = 1,NBLK
        MATDIM = MATDIM + NDIM(IBLK)*NDIM(IBLK)
    200 END DO

    AMAX = ZERO
    DO 205 I=1,MATDIM
        AMAX = MAX( AMAX, ABS(AMATRX(I)) )
    205 END DO
    IF (AMAX == ZERO) THEN
        WRITE (LVPRI,3000) NBLK
        GO TO 800
    END IF
    IF (FFMIN <= AMAX .AND. AMAX <= FFMAX) THEN
    !        use F output format
        PFMT = '(A1,I7,2X,8F15.8)'
    ELSE
    !        use 1PD output format
        PFMT = '(A1,I7,2X,1P,8D15.6)'
    END IF

    IOFF = 0
    DO 100 IBLK = 1,NBLK
        IDIM = NDIM(IBLK)
        IF (IDIM == 0) GO TO 100
        I2DIM = IDIM*IDIM
        DO 10 I=1,I2DIM
            IF (AMATRX(IOFF+I) /= ZERO) GO TO 15
        10 END DO
    WRITE (LVPRI,1100) IBLK
    GO TO 90
    15 CONTINUE
    WRITE (LVPRI,1200) IBLK
    LAST = MIN(IDIM,KCOL)
    DO 2 BEGIN = 1,IDIM,KCOL
        WRITE (LVPRI,1000) (COLUMN,I,I = BEGIN,LAST)
        DO 1 K = 1,IDIM
            IOFFK = IOFF + K
            DO 4 I=BEGIN,LAST
                IF (AMATRX(IOFFK+(I-1)*IDIM) /= ZERO) GO TO 5
            4 END DO
        GO TO 1
        5 WRITE (LVPRI,PFMT) CTL,K,(AMATRX(IOFFK+(I-1)*IDIM), &
        I = BEGIN,LAST)
        1 END DO
    LAST = MIN(LAST+KCOL,IDIM)
    2 END DO
    90 IOFF = IOFF + I2DIM
    100 END DO

    800 CONTINUE
    RETURN
    3000 FORMAT (/5X,'All',I3,' blocks zero matrices.')
    1100 FORMAT (/5X,'*** Block',I3,' zero matrix. ***')
    1200 FORMAT (/5X,'*** Block',I3,' ***')
    1000 FORMAT (/12X,6(3X,A6,I4,2X),(3X,A6,I4))
! 2000 FORMAT (A1,'Row',I4,2X,1P,8D15.6)
! 2000 FORMAT (A1,I7,2X,1P,8D15.6)
    END SUBROUTINE OUTPTB_
!  /* Deck output */
    SUBROUTINE OUTPUT_(AMATRX,ROWLOW,ROWHI,COLLOW,COLHI,ROWDIM,COLDIM, &
    NCTL,LVPRI)
!.......................................................................
! Revised 15-Dec-1983 by Hans Jorgen Aa. Jensen.
!         16-Jun-1986 hjaaj ( removed Hollerith )

! OUTPUT_ PRINTS A REAL MATRIX IN FORMATTED FORM WITH NUMBERED ROWS
! AND COLUMNS.  THE INPUT IS AS FOLLOWS;

!        AMATRX(',').........MATRIX TO BE OUTPUT_

!        ROWLOW..............ROW NUMBER AT WHICH OUTPUT_ IS TO BEGIN

!        ROWHI...............ROW NUMBER AT WHICH OUTPUT_ IS TO END

!        COLLOW..............COLUMN NUMBER AT WHICH OUTPUT_ IS TO BEGIN

!        COLHI...............COLUMN NUMBER AT WHICH OUTPUT_ IS TO END

!        ROWDIM..............ROW DIMENSION OF AMATRX(',')

!        COLDIM..............COLUMN DIMENSION OF AMATRX(',')

!        NCTL................CARRIAGE CONTROL FLAG; 1 FOR SINGLE SPACE
!                                                   2 FOR DOUBLE SPACE
!                                                   3 FOR TRIPLE SPACE
!                            hjaaj: negative for 132 col width

! THE PARAMETERS THAT FOLLOW MATRIX ARE ALL OF TYPE INTEGER*4.  THE
! PROGRAM IS SET UP TO HANDLE 5 COLUMNS/PAGE WITH A 1P,5D24.15 FORMAT
! FOR THE COLUMNS.  IF A DIFFERENT NUMBER OF COLUMNS IS REQUIRED,
! CHANGE FORMATS 1000 AND 2000, AND INITIALIZE KCOL WITH THE NEW NUMBER
! OF COLUMNS.

! AUTHOR;  NELSON H.F. BEEBE, QUANTUM THEORY PROJECT, UNIVERSITY OF
!          FLORIDA, GAINESVILLE
! REVISED; FEBRUARY 26, 1971

!.......................................................................

#include <pcm_implicit.h>
    INTEGER ::   ROWLOW,ROWHI,COLLOW,COLHI,ROWDIM,COLDIM,BEGIN,KCOL
    DIMENSION AMATRX(ROWDIM,COLDIM)
    CHARACTER(1) :: ASA(3), BLANK, CTL
    CHARACTER   PFMT*20, COLUMN*8
    PARAMETER (ZERO=0.0D00, KCOLP=5, KCOLN=8)
    PARAMETER (FFMIN=1.0D-03, FFMAX = 1.0D03)
    DATA COLUMN/'Column  '/, BLANK/' '/, ASA/' ', '0', '-'/

    IF (ROWHI < ROWLOW) GO TO 3
    IF (COLHI < COLLOW) GO TO 3

    AMAX = ZERO
    DO 10 J = COLLOW,COLHI
        DO 10 I = ROWLOW,ROWHI
            AMAX = MAX( AMAX, ABS(AMATRX(I,J)) )
    10 END DO
    IF (AMAX == ZERO) THEN
        WRITE (LVPRI,'(/T6,A)') 'Zero matrix.'
        GO TO 3
    END IF
    IF (FFMIN <= AMAX .AND. AMAX < FFMAX) THEN
    !        use F output format
        PFMT = '(A1,I7,2X,8F14.8)'
        thrpri = 0.5D-8
    ELSE
    !        use 1PD output format
        PFMT = '(A1,I7,2X,1P,8D14.6)'
        thrpri = 1.0D-8*AMAX
    END IF

    IF (NCTL < 0) THEN
        KCOL = KCOLN
    ELSE
        KCOL = KCOLP
    END IF
    MCTL = ABS(NCTL)
    IF ((MCTL <= 3) .AND. (MCTL > 0)) THEN
        CTL = ASA(MCTL)
    ELSE
        CTL = BLANK
    END IF

    LAST = MIN(COLHI,COLLOW+KCOL-1)
    DO 2 BEGIN = COLLOW,COLHI,KCOL
        WRITE (LVPRI,1000) (COLUMN,I,I = BEGIN,LAST)
        DO 1 K = ROWLOW,ROWHI
            DO 4 I = BEGIN,LAST
                IF (abs(AMATRX(K,I)) > thrpri) GO TO 5
            4 END DO
        GO TO 1
        5 WRITE (LVPRI,PFMT) CTL,K,(AMATRX(K,I), I = BEGIN,LAST)
        1 END DO
    LAST = MIN(LAST+KCOL,COLHI)
    2 END DO
    3 RETURN
    1000 FORMAT (/10X,8(4X,A6,I4))
! 2000 FORMAT (A1,'Row',I4,2X,1P,8D14.6)
! 2000 FORMAT (A1,I7,2X,1P,8D14.6)
    END SUBROUTINE OUTPUT_
    SUBROUTINE OUTPUT_x(AMATRX,ROWLOW,ROWHI,COLLOW,COLHI,ROWDIM, &
    COLDIM,NCTL,LVPRI)
!.......................................................................
! Revised 15-Dec-1983 by Hans Jorgen Aa. Jensen.
!         16-Jun-1986 hjaaj ( removed Hollerith )
! OUTPUT_x based on OUTPUT 2007, Hans Jorgen Aa. Jensen
!         prints 7F10.3/D10.3 instead of 5F14.6/D14.6

! OUTPUT_ PRINTS A REAL MATRIX IN FORMATTED FORM WITH NUMBERED ROWS
! AND COLUMNS.  THE INPUT IS AS FOLLOWS;

!        AMATRX(',').........MATRIX TO BE OUTPUT_

!        ROWLOW..............ROW NUMBER AT WHICH OUTPUT_ IS TO BEGIN

!        ROWHI...............ROW NUMBER AT WHICH OUTPUT_ IS TO END

!        COLLOW..............COLUMN NUMBER AT WHICH OUTPUT_ IS TO BEGIN

!        COLHI...............COLUMN NUMBER AT WHICH OUTPUT_ IS TO END

!        ROWDIM..............ROW DIMENSION OF AMATRX(',')

!        COLDIM..............COLUMN DIMENSION OF AMATRX(',')

!        NCTL................CARRIAGE CONTROL FLAG; 1 FOR SINGLE SPACE
!                                                   2 FOR DOUBLE SPACE
!                                                   3 FOR TRIPLE SPACE
!                            hjaaj: negative for 132 col width

! THE PARAMETERS THAT FOLLOW MATRIX ARE ALL OF TYPE INTEGER*4.  THE
! PROGRAM IS SET UP TO HANDLE 5 COLUMNS/PAGE WITH A 1P,5D24.15 FORMAT
! FOR THE COLUMNS.  IF A DIFFERENT NUMBER OF COLUMNS IS REQUIRED,
! CHANGE FORMATS 1000 AND 2000, AND INITIALIZE KCOL WITH THE NEW NUMBER
! OF COLUMNS.

! AUTHOR;  NELSON H.F. BEEBE, QUANTUM THEORY PROJECT, UNIVERSITY OF
!          FLORIDA, GAINESVILLE
! REVISED; FEBRUARY 26, 1971

!.......................................................................

#include <pcm_implicit.h>
    INTEGER ::   ROWLOW,ROWHI,COLLOW,COLHI,ROWDIM,COLDIM,BEGIN,KCOL
    DIMENSION AMATRX(ROWDIM,COLDIM)
    CHARACTER(1) :: ASA(3), BLANK, CTL
    CHARACTER   PFMT*20, COLUMN*8
    PARAMETER (ZERO=0.0D00, KCOLP=7, KCOLN=10)
    PARAMETER (FFMIN=1.0D-03, FFMAX = 1.0D03)
    DATA COLUMN/'Column  '/, BLANK/' '/, ASA/' ', '0', '-'/

    IF (ROWHI < ROWLOW) GO TO 3
    IF (COLHI < COLLOW) GO TO 3

    AMAX = ZERO
    DO 10 J = COLLOW,COLHI
        DO 10 I = ROWLOW,ROWHI
            AMAX = MAX( AMAX, ABS(AMATRX(I,J)) )
    10 END DO
    IF (AMAX == ZERO) THEN
        WRITE (LVPRI,'(/T6,A)') 'Zero matrix.'
        GO TO 3
    END IF
    IF (FFMIN <= AMAX .AND. AMAX <= FFMAX) THEN
    !        use F output format
        PFMT = '(A1,I5,2X,10F10.3)'
    ELSE
    !        use 1PD output format
        PFMT = '(A1,I5,2X,1P,10D10.3)'
    END IF

    IF (NCTL < 0) THEN
        KCOL = KCOLN
    ELSE
        KCOL = KCOLP
    END IF
    MCTL = ABS(NCTL)
    IF ((MCTL <= 3) .AND. (MCTL > 0)) THEN
        CTL = ASA(MCTL)
    ELSE
        CTL = BLANK
    END IF

    LAST = MIN(COLHI,COLLOW+KCOL-1)
    DO 2 BEGIN = COLLOW,COLHI,KCOL
        WRITE (LVPRI,1000) (COLUMN,I,I = BEGIN,LAST)
        DO 1 K = ROWLOW,ROWHI
            DO 4 I = BEGIN,LAST
                IF (abs(AMATRX(K,I)) > 0.5D-3) GO TO 5
            4 END DO
        GO TO 1
        5 WRITE (LVPRI,PFMT) CTL,K,(AMATRX(K,I), I = BEGIN,LAST)
        1 END DO
    LAST = MIN(LAST+KCOL,COLHI)
    2 END DO
    3 RETURN
    1000 FORMAT (/9X,10(A6,I3,1X),(3X,A6,I4))
! 2000 FORMAT (A1,'Row',I4,2X,1P,8D15.6)
! 2000 FORMAT (A1,I7,2X,1P,8D15.6)
    end SUBROUTINE OUTPUT_x
!  /* Deck prtab */
    SUBROUTINE PRTAB_(NTABLE,TABLE,TEXT,LVPRI)

! 28-Dec-1989 Hans Joergen Aa. Jensen

! Purpose: print tables of text (e.g. from input parsing routines)

    CHARACTER*(*) TABLE(NTABLE), TEXT
    LTEXT = LEN(TEXT)
    LTEXT = MIN(70,LTEXT)
    WRITE (LVPRI,'(//1X,A/1X,70A1/)') TEXT(1:LTEXT),('=',I=1,LTEXT)
    DO 100 I = 1,NTABLE
        IF (INDEX(TABLE(I),'XXXX') == 0) THEN
            WRITE (LVPRI,'(T6,A)') TABLE(I)
        END IF
    100 END DO
    WRITE (LVPRI,'()')
    RETURN
    END SUBROUTINE PRTAB_
!  /* Deck prvec */
    SUBROUTINE PRVEC_(NDIM,VEC,INCVEC,THRESH,MAXLIN,LUOUT)

! 19-Aug-1989 Hans Joergen Aa. Jensen

!   print VEC(1:NDIM*INCVEC:INCVEC)

!   NDIM      : Number of elements in vector VEC
!   VEC(:)    : Vector to be printed
!   INCVEC    : Increment between each element in vector VEC
!   THRESH    : Print threshold for vector with unit norm
!               (if THRESH .lt. 0 then -THRESH is used without
!                renormalization).
!   MAXLIN    : max. lines of output with vector elements
!   LUOUT     : output unit


#include <pcm_implicit.h>
    DIMENSION VEC(*)
    PARAMETER (D0 = 0.0D0, D1 = 1.0D0)
    PARAMETER (D1LOW = D1 - 1.0D-10, D1HGH = D1 + 1.0D-10)
    LOGICAL ::   VSCALE

!     Test input

    NERR = 0
    IF (NDIM   <= 0) THEN
        WRITE (LUOUT,'(/5X,A)') &
        'No print from PRVEC_ because NDIM .le. 0'
        NERR = NERR + 1
    END IF
    IF (INCVEC <= 0) THEN
        WRITE (LUOUT,'(/5X,A)') &
        'No print from PRVEC_ because INCVEC .le. 0'
        NERR = NERR + 1
    END IF
    IF (NERR > 0) RETURN



    C2NRM = DNORM2__(NDIM,VEC,INCVEC)
    IF (THRESH <= D0) THEN
        THRES2 = -THRESH
        VSCALE = .FALSE.
    ELSE
        THRES2 =  THRESH * C2NRM
        VSCALE = (C2NRM .LT. D1LOW .OR. C2NRM .GT. D1HGH)
    END IF
    IF (VSCALE) THEN
        WRITE (LUOUT,'(/2A,1P,D12.4)') &
        ' Print of vector elements (vector scaled to unit norm)', &
        ' larger than',ABS(THRESH)
        SCALE = D1 / C2NRM
    ELSE
        WRITE (LUOUT,'(/A,1P,D12.4)') &
        ' Print of vector elements larger than',ABS(THRESH)
        SCALE = D1
    END IF

    IPR   = 0
    IZER  = 0
    C2SUM = D0

    DO 300 I = 1, NDIM
        NA = 1 + (I-1)*INCVEC
        IF (ABS(VEC(NA)) <= THRES2 .OR. IPR >= MAXLIN) THEN
            C2SUM = C2SUM + VEC(NA)*VEC(NA)
            IF (VEC(NA) == D0) IZER = IZER + 1
        ELSE
            IF (MOD(IPR,5) == 0) WRITE (LUOUT,'()')
            IPR = IPR + 1
            IF (VSCALE) THEN
                WRITE(LUOUT,50)I,VEC(NA),SCALE*VEC(NA)
            ELSE
                WRITE(LUOUT,60)I,VEC(NA)
            END IF
            50 FORMAT(3X,'Element',I10,3X,'coefficient',1P,D10.2,0P, &
            3X,'scaled to unit norm',F10.6)
            60 FORMAT(6X,'Element',I12,3X,'coefficient',F20.10)
        END IF
    300 END DO
    IF (IPR >= MAXLIN) THEN
            
    !     *** We have reached the print limit
            
        WRITE (LUOUT,910) IPR
    END IF
    910 FORMAT(/' Print limit of',I6,' elements has been reached.')
    C2SUM = SQRT(C2SUM)
    IF (IZER == NDIM) THEN
        WRITE (LUOUT,920) NDIM
    ELSE
        WRITE (LUOUT,930) NDIM,NDIM-IPR,IZER,C2SUM,C2NRM
    END IF
    920 FORMAT(/' Length of vector                      :',I10, &
    /' All elements are zero.',/)
    930 FORMAT(/' Length of vector                      :',I10, &
    /' Number of elements not printed        :',I10, &
    /' Number of zero elements               :',I10, &
    /' Total norm of coefficients not printed:',F10.6, &
    /' (the coefficients are normalized to    ',F10.6,')',/)
    RETURN

! End of PRVEC_

    END SUBROUTINE PRVEC_
!  /* Deck prmgn */
    SUBROUTINE PRMGN_(NDIM,VEC,INCVEC,NPOT,LUOUT)

! 19-Aug-1989 hjaaj

!   NDIM      : Number of elements in vector VEC
!   VEC(:)    : Vector to be printed
!   INCVEC    : Increment between each element in vector VEC
!   NPOT      : IBASE**-NPOT is lowest magnitude considered
!   LUOUT     : output unit


#include <pcm_implicit.h>
    DIMENSION VEC(*)
    PARAMETER (D0 = 0.0D0, D1 = 1.0D0)
    PARAMETER (IBASE = 10)
!     IBASE     : Base for magnitude

    FACMGN = IBASE
    FACMGN = D1 / FACMGN

!     Test input

    NERR = 0
    IF (NDIM   <= 0) THEN
        WRITE (LUOUT,*) 'No print from PRMGN_ because NDIM .le. 0'
        NERR = NERR + 1
    END IF
    IF (INCVEC <= 0) THEN
        WRITE (LUOUT,*) 'No print from PRMGN_ because INCVEC .le. 0'
        NERR = NERR + 1
    END IF
    IF (NERR > 0) RETURN



    C2NRM = DNORM2__(NDIM,VEC,INCVEC)
    IZER  = 0
    NLAST = NDIM*INCVEC
    DO 100 NA = 1, NLAST, INCVEC
        IF (VEC(NA) == D0) IZER = IZER + 1
    100 END DO

!.... SIZE OF COEFFICIENTS

    WRITE(LUOUT,'(/A)')  '  Magnitude of coefficients '
    WRITE(LUOUT,'( A)')  ' ==========================='
    IF (IZER == NDIM) THEN
        WRITE(LUOUT,'(/A)') ' All elements are zero.'
        GO TO 9999
    END IF

    IF (C2NRM > (D1 + 1.0D-10) ) THEN
        WRITE (LUOUT,'(/A,1P,E23.12,A,/A)') &
        ' (Vector norm is',C2NRM,')', &
        ' (Range is scaled to a vector norm of 1)'
    END IF
    XMIN  = MAX(C2NRM,D1)
!     ... max possible coefficient is C2NRM

    WRITE (LUOUT,'(/2A,/2A)') &
    '     Range        # of elements', &
    '     Norm squared         Accumulated', &
    ' -------------   --------------', &
    '   ---------------     ---------------'

! Output will look like:

!     Range        # of elements     Norm squared         Accumulated
! -------------   --------------   ---------------     ---------------
! 10- 1 to  1                2      1.12345678E-01      1.12345678E-01

    C2NRM = D0
    IDET  = IZER
!. LOOP OVER INTER_VALS
    IPOT  = 0
    200 CONTINUE
    CLNORM = D0
    INRANG = 0
    XMAX   = XMIN
    IF (IPOT == 0) XMAX = 2*XMAX
!        ... to make certain we don't omit any because of round-off
    XMIN   = XMIN * FACMGN
    DO 300 NA = 1, NLAST, INCVEC
        IF( ABS(VEC(NA)) <= XMAX  .AND. &
        ABS(VEC(NA)) > XMIN ) THEN
            INRANG = INRANG + 1
            CLNORM = CLNORM + VEC(NA) ** 2
        END IF
    300 END DO
    C2NRM = C2NRM + CLNORM


    IF (INRANG > 0) THEN
        IF (IPOT == 0) THEN
            WRITE(LUOUT,'(I3,A,I14,1P,2E20.8)') &
            IBASE,'- 1 to  1   ',INRANG,CLNORM,C2NRM
        ELSE
            WRITE(LUOUT,'(I3,A,I2,A,I3,A,I2,I14,1P,2E20.8)') &
            IBASE,'-',IPOT+1,' to',IBASE,'-',IPOT, &
            INRANG,CLNORM,C2NRM
        END IF
    END IF


    IDET = IDET + INRANG
    IPOT = IPOT + 1
    IF( IDET < NDIM .AND. IPOT < NPOT ) GOTO 200

    ISML = NDIM - IDET
    IF (ISML > 0) WRITE (LUOUT,'(A,I13)') ' Other non-zero ',ISML
    IF (IZER > 0) WRITE (LUOUT,'(A,I13)') ' Exact zero     ',IZER
    WRITE (LUOUT,'(/A,I13/)') ' Total          ',NDIM


    9999 CONTINUE
    RETURN

! End of PRMGN_

    END SUBROUTINE PRMGN_
!  /* Deck pnzvec */
    FUNCTION PNZVEC_(NDIM,VEC,INCVEC,THRZER)

! Copyright 24-Mar-1993 Hans Joergen Aagaard Jensen
!   return % non-zero elements of VEC(1:NDIM*INCVEC:INCVEC)

!   NDIM      : Number of elements in vector VEC
!   VEC(:)    : Vector
!   INCVEC    : Increment between each element in vector VEC
!   THRZER    : Threshold for zero elements


#include <pcm_implicit.h>
    DIMENSION VEC(*)
#include <pcm_priunit.h>
    PARAMETER (D0 = 0.0D0, D100 = 100.0D0)

!     Test input

    NERR = 0
    IF (NDIM   <= 0) THEN
        WRITE (LVPRI,'(/5X,A,I10)') &
        'No percentage from PNZVEC_ because NDIM .le. 0; NDIM =',NDIM
        NERR = NERR + 1
    END IF
    IF (INCVEC <= 0) THEN
        WRITE (LVPRI,'(/5X,A,I10)') &
        'No percentage from PNZVEC_ because INCVEC .le. 0;INCVEC =', &
        INCVEC
        NERR = NERR + 1
    END IF
    IF (NERR > 0) THEN
        PNZVEC_ = -D100
        RETURN
    END IF

    THRESH = MAX(D0,THRZER)
    NNZ = 0
    DO 100 I = 1,NDIM*INCVEC,INCVEC
        IF (ABS(VEC(I)) > THRESH) NNZ = NNZ + 1
    100 END DO
    DNZ    = NNZ
    DNDIM  = NDIM
    PNZVEC_ = DNZ*D100 / DNDIM
    RETURN
    END FUNCTION PNZVEC_
!  /* Deck coutput */
    SUBROUTINE COUTPUT__(AMATRX,ROWLOW,ROWHI,COLLOW,COLHI,ROWDIM, &
    COLDIM,NCTL,LVPRI)
!.......................................................................
! Revised 15-Dec-1983 by Hans Jorgen Aa. Jensen.
!         16-Jun-1986 hjaaj ( removed Hollerith )
!         20-Mar-2001 ov    complex matrix verison

! OUTPUT_ PRINTS A REAL MATRIX IN FORMATTED FORM WITH NUMBERED ROWS
! AND COLUMNS.  THE INPUT IS AS FOLLOWS;

!        AMATRX(',').........MATRIX TO BE OUTPUT_

!        ROWLOW..............ROW NUMBER AT WHICH OUTPUT_ IS TO BEGIN

!        ROWHI...............ROW NUMBER AT WHICH OUTPUT_ IS TO END

!        COLLOW..............COLUMN NUMBER AT WHICH OUTPUT_ IS TO BEGIN

!        COLHI...............COLUMN NUMBER AT WHICH OUTPUT_ IS TO END

!        ROWDIM..............ROW DIMENSION OF AMATRX(',')

!        COLDIM..............COLUMN DIMENSION OF AMATRX(',')

!        NCTL................CARRIAGE CONTROL FLAG; 1 FOR SINGLE SPACE
!                                                   2 FOR DOUBLE SPACE
!                                                   3 FOR TRIPLE SPACE

! THE PARAMETERS THAT FOLLOW MATRIX ARE ALL OF TYPE INTEGER*4.  THE
! PROGRAM IS SET UP TO HANDLE 5 COLUMNS/PAGE WITH A 1P,5D24.15 FORMAT
! FOR THE COLUMNS.  IF A DIFFERENT NUMBER OF COLUMNS IS REQUIRED,
! CHANGE FORMATS 1000 AND 2000, AND INITIALIZE KCOL WITH THE NEW NUMBER
! OF COLUMNS.

! AUTHOR;  NELSON H.F. BEEBE, QUANTUM THEORY PROJECT, UNIVERSITY OF
!          FLORIDA, GAINESVILLE
! REVISED; FEBRUARY 26, 1971

!.......................................................................

    IMPLICIT DOUBLE COMPLEX (A-H,O-Z)
    INTEGER ::   ROWLOW,ROWHI,COLLOW,COLHI,ROWDIM,COLDIM,BEGIN,KCOL
    DIMENSION AMATRX(ROWDIM,COLDIM)
    CHARACTER(1) :: ASA(3), BLANK, CTL
    CHARACTER   PFMT*20, COLUMN*8
    DOUBLE PRECISION :: ZERO, FFMIN, FFMAX
    PARAMETER (ZERO=0.D00, KCOLP=4, KCOLN=6)
    PARAMETER (FFMIN=1.D-3, FFMAX = 1.D3)
    DATA COLUMN/'Column  '/, BLANK/' '/, ASA/' ', '0', '-'/
    REAL(8) :: AMAX

    IF (ROWHI < ROWLOW) GO TO 3
    IF (COLHI < COLLOW) GO TO 3

    AMAX = ZERO
    DO 10 J = COLLOW,COLHI
        DO 10 I = ROWLOW,ROWHI
            AMAX = MAX( AMAX, ABS(AMATRX(I,J)) )
    10 END DO
    IF (AMAX == ZERO) THEN
        WRITE (LVPRI,'(/T6,A)') 'Zero matrix.'
        GO TO 3
    END IF
    IF (FFMIN <= AMAX .AND. AMAX <= FFMAX) THEN
    !        use F output format
        PFMT = '(A1,I3,2X,4(2F10.6))'
    ELSE
    !        use 1PD output format
        PFMT = '(A1,I7,2X,1P,4(2D12.6))'
    END IF

    IF (NCTL < 0) THEN
        KCOL = KCOLN
    ELSE
        KCOL = KCOLP
    END IF
    MCTL = ABS(NCTL)
    IF ((MCTL <= 3) .AND. (MCTL > 0)) THEN
        CTL = ASA(MCTL)
    ELSE
        CTL = BLANK
    END IF

    LAST = MIN(COLHI,COLLOW+KCOL-1)
    DO 2 BEGIN = COLLOW,COLHI,KCOL
        WRITE (LVPRI,1000) (COLUMN,I,I = BEGIN,LAST)
        DO 1 K = ROWLOW,ROWHI
            DO 4 I = BEGIN,LAST
                IF (AMATRX(K,I) /= ZERO) GO TO 5
            4 END DO
        GO TO 1
        5 WRITE (LVPRI,PFMT) CTL,K,(AMATRX(K,I), I = BEGIN,LAST)
        1 END DO
    LAST = MIN(LAST+KCOL,COLHI)
    2 END DO
    3 RETURN
    1000 FORMAT (/5X,6(3X,A6,I4,7X),(3X,A6,I4))
! 2000 FORMAT (A1,'Row',I4,2X,1P,8D15.6)
! 2000 FORMAT (A1,I7,2X,1P,8D15.6)
    END SUBROUTINE COUTPUT__
! --- end of printpkg.F ---
