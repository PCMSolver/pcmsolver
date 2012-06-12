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
      subroutine makecav_(xtscor_,ytscor_,ztscor_,ar_,
     & xsphcor_,ysphcor_, zsphcor_, rsph_, nts_, nesfp_, xe_, ye_, ze_,
     & rin_, avgarea_, rsolv_, work2, lwork2)

#include <pcm_implicit.h>
#include <pcm_maxorb.h>
#include <pcm_maxaqn.h>
#include <pcm_priunit.h>
#include <pcm_iratdef.h>
#include <pcm_pcmdef.h>
#include <pcm_mxcent.h>
#include <pcm_infpri.h>
#include <pcm_pcm.h>
#include <pcm_pcmlog.h>
#include <pcm_symmet.h>
#include <pcm_pgroup.h>

      dimension work2(*)
      double precision xtscor_(*), ytscor_(*), ztscor_(*)
      double precision xsphcor_(*), ysphcor_(*), zsphcor_(*), rsph_(*)
      double precision ar_(*), xe_(*), ye_(*), ze_(*), rin_(*)
      double precision avgarea_, rsolv_
      integer nts_, nesfp_
      logical pedra_file_exists

      LVPRI = 121201
      inquire(file = 'PEDRA.OUT', exist = pedra_file_exists)
      if (pedra_file_exists) then
         open(lvpri, 
     &        file = 'PEDRA.OUT', 
     &        status = 'unknown', 
     &        form = 'formatted',
     &        access = 'sequential')
         close(lvpri, status = 'delete')
      end if
      open(lvpri,
     &     file = 'PEDRA.OUT',
     &     status = 'new',
     &     form = 'formatted',
     &     access = 'sequential')
      rewind(lvpri)
      MAXREP=0
      NSYM=1
      AREATS=avgArea_
      ICESPH=1
      IPRPCM=3
      GROUP = 'C1'
      RSOLV = rsolv_
      RET = 100.0D0
      NESFP = nesfp_
      do i=1, NESFP
         XE(i) = xe_(i)
         YE(i) = ye_(i)
         ZE(i) = ze_(i)
         RIN(i) = rin_(i)
         ALPHA(i)=1.0d0
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
      
      CALL PEDRA_M_(WORK2, LWORK2)

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

      close(lvpri)

      return
      END
