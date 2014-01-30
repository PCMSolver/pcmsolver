    module pedra_print

    implicit none

    public output

    contains

!...   Dalton, Release 2.0, pdpack/printpkg.F
!...
!...   These routines are in the public domain and can be
!...   used freely in other programs.
!...

    SUBROUTINE OUTPUT(AMATRX,ROWLOW,ROWHI,COLLOW,COLHI,ROWDIM,COLDIM, &
    NCTL,LVPRI)
!.......................................................................
! Revised 15-Dec-1983 by Hans Jorgen Aa. Jensen.
!         16-Jun-1986 hjaaj ( removed Hollerith )

! OUTPUT PRINTS A REAL MATRIX IN FORMATTED FORM WITH NUMBERED ROWS
! AND COLUMNS.  THE INPUT IS AS FOLLOWS;

!        AMATRX(',').........MATRIX TO BE OUTPUT

!        ROWLOW..............ROW NUMBER AT WHICH OUTPUT IS TO BEGIN

!        ROWHI...............ROW NUMBER AT WHICH OUTPUT IS TO END

!        COLLOW..............COLUMN NUMBER AT WHICH OUTPUT IS TO BEGIN

!        COLHI...............COLUMN NUMBER AT WHICH OUTPUT IS TO END

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

    INTEGER, intent(in) :: ROWLOW,ROWHI,COLLOW,COLHI,ROWDIM,COLDIM
    integer, intent(in) :: nctl, lvpri
    real(8), intent(in) :: AMATRX(ROWDIM,COLDIM)
    integer :: BEGIN,KCOL
    CHARACTER(1) :: ASA(3), BLANK, CTL
    CHARACTER   PFMT*20, COLUMN*8
    real(8), parameter :: zero = 0.0d0, ffmin = 1.0d-03, ffmax = 1.0d03
    integer, parameter :: kcolp = 5, kcoln = 8
    DATA COLUMN/'Column  '/, BLANK/' '/, ASA/' ', '0', '-'/

    real(8) :: amax, thrpri
    integer :: i, j, k, mctl, last

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
    END SUBROUTINE OUTPUT

    end module
