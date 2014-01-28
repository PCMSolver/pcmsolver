    SUBROUTINE DZERO_(DX,LENGTH)
#include "pcm_implicit.h"

! Last revision 5-May-1984 by Hans Jorgen Aa. Jensen

!   Subroutine DZERO_ sets a real array of length *LENGTH*
!   to zero.
!...................................................................
    DIMENSION DX(*)

    IF (LENGTH <= 0) RETURN

    DO I = 1,LENGTH
        DX(I) = 0.0D00
    END DO

    RETURN
    END SUBROUTINE DZERO_

    SUBROUTINE AROUND_(HEAD)
    CHARACTER HEAD*(*)
#include <pcm_priunit.h>
    LHEAD  = LNBLNK(HEAD)
    LNG    = LHEAD + 2
    IND = MAX(1,(80 - LNG)/2 + 1)
    WRITE (LVPRI,'(//150A)') (' ',I=1,IND), '+', ('-',I=1,LNG), '+'
    WRITE (LVPRI,'(150A)')   (' ',I=1,IND), '! ', HEAD(1:LHEAD), ' !'
    WRITE (LVPRI,'(150A)')   (' ',I=1,IND), '+', ('-',I=1,LNG), '+'
!x    WRITE (LVPRI,'(//150A)') (' ',I=1,IND), '.', ('-',I=1,LNG), '.'
!x    WRITE (LVPRI,'(150A)')   (' ',I=1,IND), '| ', HEAD(1:LHEAD), ' |'
!x    WRITE (LVPRI,'(150A)')   (' ',I=1,IND), '`', ('-',I=1,LNG), ''''
    WRITE (LVPRI,'()')
    RETURN
    END SUBROUTINE AROUND_

    FUNCTION DNORM2__(N,DX,INCX)

!     Forms the two-norm of a vector.
! 19-Sep-1988 -- hjaaj -- based on DNRM2_ from LINPACK
!     This version does not use extended precision for intermediates
!     as the LINPACK version does.C 1) DNORM2__ (emulate ESSL DNORM2: do not use \
! xtended precision for intermediates)

!     Equivalent to DNORM2__ in IBM's ESSL library.

!     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.
!     DNRM2_: JACK DONGARRA, LINPACK, 3/11/78.

#include <pcm_implicit.h>

    DIMENSION DX(*)
    PARAMETER ( ZERO = 0.0D0 )

    IF (N <= 0) THEN
        DNORM2__ = ZERO
        RETURN
    END IF
    DTEMP  = ZERO
    IF(INCX == 1)GO TO 20

!        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
!          NOT EQUAL TO 1

    IX = 1
    IF(INCX < 0)IX = (-N+1)*INCX + 1
    DO 10 I = 1,N
        DTEMP = DTEMP + DX(IX)*DX(IX)
        IX = IX + INCX
    10 END DO
    DNORM2__ = SQRT(DTEMP)
    RETURN

!        CODE FOR BOTH INCREMENTS EQUAL TO 1


!        CLEAN-UP LOOP
!                                       20 M = MOD(N,5)
    20 M = MOD(N,5)
    IF( M == 0 ) GO TO 40
    DO 30 I = 1,M
        DTEMP = DTEMP + DX(I)*DX(I)
    30 END DO
    IF( N < 5 ) GO TO 60
    40 MP1 = M + 1
    DO 50 I = MP1,N,5
        DTEMP = DTEMP + DX(I)*DX(I) + DX(I + 1)*DX(I + 1) + &
        DX(I + 2)*DX(I + 2) + DX(I + 3)*DX(I + 3) + DX(I + 4)*DX(I + 4)
    50 END DO
    60 DNORM2__ = SQRT(DTEMP)
    RETURN
    END FUNCTION DNORM2__

    SUBROUTINE ERRWRK_ (STRING,LNEED,LAVAIL)

! Version 6-Mar-1985 by hjaaj

    CHARACTER*(*) STRING

#include <pcm_priunit.h>

! m      CALL QENTER_('ERRWRK_')
    IF (LNEED >= 0) THEN
        WRITE (LVPRI,1010) STRING,LNEED,LAVAIL
    ELSE
        WRITE (LVPRI,1020) STRING,-LNEED,LAVAIL
    END IF
! m      CALL QTRACE_(LVPRI)
    STOP

    1010 FORMAT(/'  FATAL ERROR, insufficient core for ',A, &
    /T16,'( Need:',I10,', available:',I10,' )')
    1020 FORMAT(/'  FATAL ERROR, insufficient core for ',A, &
    //T16,'Need      :',I10,' + uncalculated amount', &
    /T16,'Available :',I10)
    END SUBROUTINE ERRWRK_

    subroutine get_point_group(int_pgroup, char_pgroup)
        
        integer          :: int_pgroup
        character(len=3) :: char_pgroup
 
        if (int_group == 0) then
                char_pgroup = 'C1'
        else if (int_group == 1) then
                char_pgroup = 'Cs'
        else if (int_group == 2) then
                char_pgroup = 'Ci'
        else if (int_group == 3) then
                char_pgroup = 'C2'
        else if (int_group == 4) then
                char_pgroup = 'D2'
        else if (int_group == 5) then
                char_pgroup = 'C2v'
        else
                char_pgroup = 'D2h'
        end if

    end subroutine get_point_group
