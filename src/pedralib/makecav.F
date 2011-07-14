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
     & rin_, avgarea_)

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
      double precision avgarea_
      integer nts_, nesfp_


      LUPRI = 3
      MAXREP=0
      NSYM=1
      AREATS=avgArea_
      ICESPH=1
      IPRPCM=20
      GROUP = 'C1'
      RSOLV = 0.0D0
      RET = 100.0D0
      print *,"avgearea",  avgarea_
      NESFP = nesfp_
      do i=1, NESFP
         XE(i) = xe_(i)
         YE(i) = ye_(i)
         ZE(i) = ze_(i)
         RIN(i) = rin_(i)
         ALPHA(i)=1
      enddo
      
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
      
      nts_=nts
      do i=1,NTS
         xtscor_(i) = xtscor(i)
         ytscor_(i) = ytscor(i)
         ztscor_(i) = ztscor(i)
         ar_(i) = as(i)
         xsphcor_(i) = xe(isphe(i))
         ysphcor_(i) = ye(isphe(i))
         zsphcor_(i) = ze(isphe(i))
         rsph_(i) = re(isphe(i))
      enddo
      return
      END
