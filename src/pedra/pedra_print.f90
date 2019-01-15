!
! PCMSolver, an API for the Polarizable Continuum Model
! Copyright (C) 2019 Roberto Di Remigio, Luca Frediani and contributors.
!
! This file is part of PCMSolver.
!
! PCMSolver is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! PCMSolver is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with PCMSolver.  If not, see <http://www.gnu.org/licenses/>.
!
! For information on the complete list of contributors to the
! PCMSolver API, see: <http://pcmsolver.readthedocs.io/>
!

module pedra_print

use pedra_precision

implicit none

public output

private

contains

    ! Dalton, Release 2.0, pdpack/printpkg.F
    !
    ! These routines are in the public domain and can be
    ! used freely in other programs.
    ! Copied to PEDRA by Roberto Di Remigio, 2013
    subroutine output(amatrx, rowlow, rowhi, collow, colhi, rowdim, coldim, nctl, lvpri)
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
    real(kind=dp),          intent(in) :: amatrx(rowdim, coldim)
    integer(kind=regint_k) :: begin, kcol
    character(1) :: asa(3), blank, ctl
    character   pfmt*20, column*8
    real(kind=dp), parameter :: zero = 0.0_dp, ffmin = 1.0e-3_dp, ffmax = 1.0e3_dp
    integer(kind=regint_k), parameter :: kcolp = 5, kcoln = 8
    DATA COLUMN/'Column  '/, BLANK/' '/, ASA/' ', '0', '-'/

    real(kind=dp) :: amax, thrpri
    integer(kind=regint_k) :: i, j, k, l, mctl, last

    if (rowhi < rowlow) return
    if (colhi < collow) return

    amax = zero
    do j = collow,colhi
        do i = rowlow,rowhi
            amax = max(amax, abs(amatrx(i, j)))
        end do
    end do
    if (amax == zero) then
        write (lvpri,'(/t6,a)') 'Zero matrix'
        return
    end if
    if (ffmin <= amax .and. amax < ffmax) then
    !        use F output format
        pfmt = '(a1,i7,2x,8f14.8)'
        thrpri = 0.5e-8_dp
    else
    !        use 1PD output format
        pfmt = '(a1,i7,2x,1p,8d14.6)'
        thrpri = 1.0e-8_dp * amax
    end if

    if (nctl < 0) then
        kcol = kcoln
    else
        kcol = kcolp
    end if
    mctl = abs(nctl)
    if ((mctl <= 3) .and. (mctl > 0)) then
        ctl = asa(mctl)
    else
        ctl = blank
    end if

    last = min(colhi, collow + kcol - 1_regint_k)
    do begin = collow,colhi,kcol
        write(lvpri, 1000) (column, i, i = begin, last)
        do k = rowlow, rowhi
            do i = begin, last
                if (abs(amatrx(k, i)) > thrpri) then
                     write (lvpri, pfmt) ctl, k, (amatrx(k, l), l = begin, last)
                     exit
                end if
            end do
        end do
    last = min(last + kcol, colhi)
    end do
    1000 format (/10x,8(4x,a6,i4))
    end subroutine output

end module
