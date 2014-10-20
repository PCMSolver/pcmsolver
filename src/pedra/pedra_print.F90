!pcmsolver_copyright_start
!       PCMSolver, an API for the Polarizable Continuum Model
!       Copyright (C) 2013 Roberto Di Remigio, Luca Frediani and contributors
!       
!       This file is part of PCMSolver.
! 
!       PCMSolver is free software: you can redistribute it and/or modify       
!       it under the terms of the GNU Lesser General Public License as published by
!       the Free Software Foundation, either version 3 of the License, or
!       (at your option) any later version.
!                                                                            
!       PCMSolver is distributed in the hope that it will be useful,
!       but WITHOUT ANY WARRANTY; without even the implied warranty of
!       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!       GNU Lesser General Public License for more details.
!                                                                            
!       You should have received a copy of the GNU Lesser General Public License
!       along with PCMSolver.  If not, see <http://www.gnu.org/licenses/>.
! 
!       For information on the complete list of contributors to the
!       PCMSolver API, see: <http://pcmsolver.github.io/pcmsolver-doc>
!pcmsolver_copyright_end

    module pedra_print
    
    use pedra_precision

    implicit none

    public output

    private

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

    integer(kind=regint_k), intent(in) :: rowlow, rowhi 
    integer(kind=regint_k), intent(in) :: collow, colhi
    integer(kind=regint_k), intent(in) :: rowdim, coldim
    integer(kind=regint_k), intent(in) :: nctl, lvpri
    real(kind=dp),                intent(in) :: amatrx(rowdim, coldim)
    integer(kind=regint_k) :: begin, kcol
    CHARACTER(1) :: ASA(3), BLANK, CTL
    CHARACTER   PFMT*20, COLUMN*8
    real(kind=dp), parameter :: zero = 0.0d0, ffmin = 1.0d-03, ffmax = 1.0d03
    integer(kind=regint_k), parameter :: kcolp = 5, kcoln = 8
    DATA COLUMN/'Column  '/, BLANK/' '/, ASA/' ', '0', '-'/

    real(kind=dp) :: amax, thrpri
    integer(kind=regint_k) :: i, j, k, mctl, last

    IF (ROWHI < ROWLOW) go to 3
    IF (COLHI < COLLOW) go to 3

    AMAX = ZERO
    DO J = COLLOW,COLHI
        DO I = ROWLOW,ROWHI
            AMAX = MAX( AMAX, ABS(AMATRX(I,J)) )
        end do
    end do
    IF (AMAX == ZERO) THEN
        WRITE (LVPRI,'(/T6,A)') 'Zero matrix.'
        go to 3
    end if
    IF (FFMIN <= AMAX .AND. AMAX < FFMAX) THEN
    !        use F output format
        PFMT = '(A1,I7,2X,8F14.8)'
        thrpri = 0.5D-8
    ELSE
    !        use 1PD output format
        PFMT = '(A1,I7,2X,1P,8D14.6)'
        thrpri = 1.0D-8*AMAX
    end if

    IF (NCTL < 0) THEN
        KCOL = KCOLN
    ELSE
        KCOL = KCOLP
    end if
    MCTL = ABS(NCTL)
    IF ((MCTL <= 3) .AND. (MCTL > 0)) THEN
        CTL = ASA(MCTL)
    ELSE
        CTL = BLANK
    end if

    LAST = MIN(COLHI,COLLOW+KCOL-1_regint_k)
    DO 2 BEGIN = COLLOW,COLHI,KCOL
        WRITE (LVPRI,1000) (COLUMN,I,I = BEGIN,LAST)
        DO 1 K = ROWLOW,ROWHI
            DO 4 I = BEGIN,LAST
                IF (abs(AMATRX(K,I)) > thrpri) go to 5
            4 end do
        go to 1
        5 WRITE (LVPRI,PFMT) CTL,K,(AMATRX(K,I), I = BEGIN,LAST)
        1 end do
    LAST = MIN(LAST+KCOL,COLHI)
    2 end do
    3 RETURN
    1000 FORMAT (/10X,8(4X,A6,I4))
! 2000 FORMAT (A1,'Row',I4,2X,1P,8D14.6)
! 2000 FORMAT (A1,I7,2X,1P,8D14.6)
    END SUBROUTINE OUTPUT

    end module
