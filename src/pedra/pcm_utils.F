      SUBROUTINE DZERO_(DX,LENGTH)
#include "pcm_implicit.h"
C                                                                               
C Last revision 5-May-1984 by Hans Jorgen Aa. Jensen                            
C                                                                               
C   Subroutine DZERO_ sets a real array of length *LENGTH*                       
C   to zero.                                                                    
C...................................................................            
      DIMENSION DX(*)
C                                                                               
      IF (LENGTH.LE.0) RETURN
C                                                                               
      DO I = 1,LENGTH
         DX(I) = 0.0D00
      END DO
C                                                                               
      RETURN
      END

      SUBROUTINE AROUND_(HEAD)
      CHARACTER HEAD*(*)
#include <pcm_priunit.h>
      LHEAD  = LNBLNK(HEAD)
      LNG    = LHEAD + 2
      IND = MAX(1,(80 - LNG)/2 + 1)
      WRITE (LVPRI,'(//150A)') (' ',I=1,IND), '+', ('-',I=1,LNG), '+'
      WRITE (LVPRI,'(150A)')   (' ',I=1,IND), '! ', HEAD(1:LHEAD), ' !'
      WRITE (LVPRI,'(150A)')   (' ',I=1,IND), '+', ('-',I=1,LNG), '+'
Cx    WRITE (LVPRI,'(//150A)') (' ',I=1,IND), '.', ('-',I=1,LNG), '.'           
Cx    WRITE (LVPRI,'(150A)')   (' ',I=1,IND), '| ', HEAD(1:LHEAD), ' |'         
Cx    WRITE (LVPRI,'(150A)')   (' ',I=1,IND), '`', ('-',I=1,LNG), ''''          
      WRITE (LVPRI,'()')
      RETURN
      END

      FUNCTION DNORM2__(N,DX,INCX)
C                                                                               
C     Forms the two-norm of a vector.                                           
C 19-Sep-1988 -- hjaaj -- based on DNRM2_ from LINPACK                           
C     This version does not use extended precision for intermediates            
C     as the LINPACK version does.C 1) DNORM2__ (emulate ESSL DNORM2: do not use \
Cextended precision for intermediates)                                           

C     Equivalent to DNORM2__ in IBM's ESSL library.                               
C                                                                               
C     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.                          
C     DNRM2_: JACK DONGARRA, LINPACK, 3/11/78.                                   
C                                                                               
#include <pcm_implicit.h>
C                                                                               
      DIMENSION DX(*)
      PARAMETER ( ZERO = 0.0D0 )
C                                                                               
      IF (N.LE.0) THEN
         DNORM2__ = ZERO
         RETURN
      END IF
      DTEMP  = ZERO
      IF(INCX.EQ.1)GO TO 20
C                                                                               
C        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS                        
C          NOT EQUAL TO 1                                                       
C                                                                               
      IX = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      DO 10 I = 1,N
        DTEMP = DTEMP + DX(IX)*DX(IX)
        IX = IX + INCX
 10      CONTINUE
      DNORM2__ = SQRT(DTEMP)
      RETURN
C                                                                               
C        CODE FOR BOTH INCREMENTS EQUAL TO 1                                    
C                                                                               
C                                                                               
C        CLEAN-UP LOOP                                                          
C                                       20 M = MOD(N,5)
 20   M = MOD(N,5)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DTEMP = DTEMP + DX(I)*DX(I)
 30      CONTINUE
      IF( N .LT. 5 ) GO TO 60
 40    MP1 = M + 1
      DO 50 I = MP1,N,5
        DTEMP = DTEMP + DX(I)*DX(I) + DX(I + 1)*DX(I + 1) +
     *   DX(I + 2)*DX(I + 2) + DX(I + 3)*DX(I + 3) + DX(I + 4)*DX(I + 4)
 50      CONTINUE
 60       DNORM2__ = SQRT(DTEMP)
      RETURN
      END

      SUBROUTINE ERRWRK_ (STRING,LNEED,LAVAIL)
C                                                                               
C Version 6-Mar-1985 by hjaaj                                                   
C                                                                               
      CHARACTER*(*) STRING
C                                                                               
#include <pcm_priunit.h>
C                                                                               
Ckm      CALL QENTER_('ERRWRK_')
      IF (LNEED .GE. 0) THEN
         WRITE (LVPRI,1010) STRING,LNEED,LAVAIL
      ELSE
         WRITE (LVPRI,1020) STRING,-LNEED,LAVAIL
      END IF
Ckm      CALL QTRACE_(LVPRI)
      STOP
C                                                                               
 1010  FORMAT(/'  FATAL ERROR, insufficient core for ',A,
     *       /T16,'( Need:',I10,', available:',I10,' )')
 1020   FORMAT(/'  FATAL ERROR, insufficient core for ',A,
     *      //T16,'Need      :',I10,' + uncalculated amount',
     *       /T16,'Available :',I10)
      END
