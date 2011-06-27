C
C Simple input/output file for cavity class.
C (badly edited for pcm wrapper)
C
C written by Krzysztof Mozgawa, 2011
C
C Clog lines provide rather unfriendly loging feature for data loading when Creader is used
C
C Cwriter lines provide writing cavity information to file cavity.out giving data in following order:
C number of tessera (only once)
C X of tessera, Y of Tessera, Z of Tessera, area of tessera, x of corresponding sphere, 
C y of corresponding sphere, z of corresponding sphere, radius of coressponding spheres.
C     Cwriter works regardless of Clog and Creader.
C
C Creader lines provide cavity input reading from the file cavity.inp.
C     Creader works regardless of Clog and Creader.
C
C if you don't want to use external (c++) input comment out lines between Cex
C
      subroutine makecav_(work,lwork,xtscor_,ytscor_,ztscor_,ar_,
     & xsphcor_,ysphcor_, zsphcor_, rsph_, nts_, nesfp_, xe_, ye_, ze_,
     & rin_)

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

      DIMENSION WORK(*)

      double precision xtscor_(*), ytscor_(*), ztscor_(*)
      double precision xsphcor_(*), ysphcor_(*), zsphcor_(*), rsph_(*)
      double precision ar_(*), xe_(*), ye_(*), ze_(*), rin_(*)
      double precision nts_, nesfp_


Creader      open(1,file='cavity.inp', STATUS='OLD', IOSTAT=INFO)
      LUPRI = 3
Clog      open(LUPRI,file='cavity.log', STATUS='REPLACE')
      MAXREP=0
      NSYM=1
      AREATS=100000000.05
      ICESPH=1
      IPRPCM=20
      GROUP = 'C1'
      RSOLV = 0.0D0
      RET = 100.0D0
C      NESF=2

C
C     In principle this should exit the program if there is no cavity  input file, but it's not working. Therefore, it's commented out.
C
Clog      IF (INFO.EQ.1) then
Clog        write(LUPRI,*) 'Input cavity file does not exist. Exiting...'
clog        STOP
clog      end if
Creader      read(1,*) NESFP

Cex
      NESFP = nesfp_
Cex

Clog      write(LUPRI,*) 'Number of spheres=', NESFP
      do i=1, NESFP
Clog         write(LUPRI,*) 'Reading coordinates or sphere nr', i
Creader         read(1,*) XE(i), YE(i), ZE(i), RIN(i)

Cex
         XE(i) = xe_(i)
         YE(i) = ye_(i)
         ZE(i) = ze_(i)
         RIN(i) = rin_(i)
Cex

         ALPHA(i)=1
Clog         write(LUPRI,*) 'x=',XE(i),'y=',YE(i), 'z=', ZE(i),'rad=',RIN(I)
      enddo
Clog      write(LUPRI,*) 'Succesfully read the file !'
Creader      close(1)
      
      
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
      
Cwriter      open(2, file='cavity.out', STATUS="REPLACE", Iostat=INFO)
      nts_=nts
Cwriter      write(2,*) NTS
      do i=1,NTS
Cwriter         write(2,*) XTSCOR(i), YTSCOR(i), ZTSCOR(i), AS(i), 
Cweiter     &        XE(ISPHE(i)), YE(ISPHE(i)), ZE(ISPHE(i)), RE(ISPHE(i))

         xtscor_(i) = xtscor(i)
         ytscor_(i) = ytscor(i)
         ztscor_(i) = ztscor(i)
         ar_(i) = as(i)
         xsphcor_(i) = xe(isphe(i))
         ysphcor_(i) = ye(isphe(i))
         zsphcor_(i) = ze(isphe(i))
         rsph_(i) = re(isphe(i))
      enddo
Cwriter      close(2)

Clog      close(LUPRI)
      return
      END
