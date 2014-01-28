!
!     simple input reader for cavity generator
!     written by Krzysztof Mozgawa, 2010
!      
!     RDR, 280114. Put things in makecav.F inside here directly.
!
      subroutine generatecavity_cpp(xtscor_, ytscor_, ztscor_, ar_, xsphcor_,     &
          ysphcor_, zsphcor_, rsph_, nts_, nesfp_, xe_, ye_, ze_, rin_, avgArea_, &
          rsolv_, work_, lwork_) bind(c, name='generatecavity_cpp')

      use, intrinsic :: iso_c_binding    

      implicit none

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

      real(c_double)  :: xtscor_(*), ytscor_(*), ztscor_(*)
      real(c_double)  :: xsphcor_(*), ysphcor_(*), zsphcor_(*), rsph_(*)
      real(c_double)  :: ar_(*), xe_(*), ye_(*), ze_(*), rin_(*)
      real(c_double)  :: avgArea_, rsolv_, work_(*)
      integer(c_int)  :: nts_, nesfp_, lwork_
      logical(c_bool) :: pedra_file_exists

      integer(c_int)  :: i, nsym
      
      lvpri = 121201
      inquire(file = 'PEDRA.OUT', exist = pedra_file_exists)
      if (pedra_file_exists) then
         open(lvpri,                   & 
             file = 'PEDRA.OUT',       &  
             status = 'unknown',       &   
             form = 'formatted',       & 
             access = 'sequential')
         close(lvpri, status = 'delete')
      end if
      open(lvpri,                      &
          file = 'PEDRA.OUT',          &
          status = 'new',              &
          form = 'formatted',          &
          access = 'sequential')
      rewind(lvpri)
      maxrep = 0
      nsym = 1
      areats = avgArea_
      icesph = 1
      iprpcm = 3
      group = 'C1'
      rsolv = rsolv_
      omega = 40.0d+00
      fro = 0.7d+00
! ret is the minimum radius of added spheres
      ret = 0.3779452249130124d0 ! This is 0.2 ang in Bohr 
      nesfp = nesfp_
      do i = 1, nesfp
         xe(i) = xe_(i)
         ye(i) = ye_(i)
         ze(i) = ze_(i)
         rin(i) = rin_(i)
         alpha(i) = 1.0d0
      enddo
      
      nesf = nesfp

      pt(0) =  1
      pt(1) = -1
      pt(2) = -1
      pt(3) =  1
      pt(4) = -1
      pt(5) =  1
      pt(6) =  1
      pt(7) =  1
      
      call pedra_m_(work_, lwork_)

      nts_ = nts
      do i = 1, nts
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

      end subroutine generatecavity_cpp
