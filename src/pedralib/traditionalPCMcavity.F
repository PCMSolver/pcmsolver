C
C     simple input reader for cavity generator
C     written by Krzysztof Mozgawa, 2010
C
      SUBROUTINE traditional_pcm_cavity(TAREA)

#include <implicit.h>
#include <maxorb.h>
#include <maxaqn.h>
#include <priunit.h>
#include <iratdef.h>
#include <pcmdef.h>
#include <mxcent.h>
#include <infpri.h>
#include <pcm.h>
#include <pcmlog.h>
#include <symmet.h>
#include <pgroup.h>

      integer INFO, i, LWORK
      PARAMETER (LWORK=10000000)
      real WORK
      real*8 TAREA
      DIMENSION WORK(LWORK)
      integer OLDLUPRI
      
      OLDLUPRI = LUPRI
      open(5700,file='cavity.inp', STATUS='OLD', IOSTAT=INFO)
      LUPRI = 5703      
      open(LUPRI,file='cavity_generator.out', STATUS='REPLACE')
      MAXREP=0
      NSYM=1
      AREATS=TAREA
      ICESPH=1
      IPRPCM=20
      GROUP = 'C1'
      RSOLV = 0.0D0
      RET = 100.0D0
C      NESF=2
      IF (INFO.EQ.1) then
        write(LUPRI,*) 'Input cavity file does not exist. Exiting...'
        STOP
      end if
      read(5700,*) NESFP
      write(LUPRI,*) 'Number of spheres=', NESFP
      do i=1, NESFP
         write(LUPRI,*) 'Reading coordinates or sphere nr', i
         read(5700,*) XE(i), YE(i), ZE(i), RIN(i)
         ALPHA(i)=1
         write(LUPRI,*) 'x=', XE(i),'y=', YE(i), 'z=',  ZE(i),
     &        'rad=', RIN(I)
      enddo
      write(LUPRI,*) 'Succesfully read the file !'
      close(5700)
      
      NESF=NESFP

      PT(0) =  1
      PT(1) = -1
      PT(2) = -1
      PT(3) =  1
      PT(4) = -1
      PT(5) =  1
      PT(6) =  1
      PT(7) =  1
         
      CALL PEDRA_M_(WORK, LWORK)
      
      open(5702, file='output_cav', STATUS="REPLACE", Iostat=INFO)
      write(5702,*) NTS
      do i=1,NTS
         write(5702,*) XTSCOR(i), YTSCOR(i), ZTSCOR(i), AS(i), 
     &        XE(ISPHE(i)), YE(ISPHE(i)), ZE(ISPHE(i)), RIN(ISPHE(i))
      enddo
      close(5702)   
      close(LUPRI)
      LUPRI = OLDLUPRI
      END
