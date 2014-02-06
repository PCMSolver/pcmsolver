!
!     simple input reader for cavity generator
!     written by Krzysztof Mozgawa, 2010
!      
!     RDR, 280114. Put things in makecav.F inside here directly.
!
      subroutine generatecavity_cpp(maxts_, maxsph_, maxvert_,         &
          xtscor_, ytscor_, ztscor_, ar_,                              &
          xsphcor_, ysphcor_, zsphcor_, rsph_,                         &
          nts_, ntsirr_, nesfp_, addsph_,                              &
          xe_, ye_, ze_, rin_, avgArea_, rsolv_, ret_, pgroup_) bind(c, name='generatecavity_cpp')

      use, intrinsic :: iso_c_binding    
      use pedra_symmetry, only: get_point_group, point_group
      use pedra_cavity, only: polyhedra_driver

      implicit none

#include "pcm_pcmdef.h"
#include "pcm_mxcent.h"
#include "pcm_pcm.h"

      integer(c_int)  :: maxts_, maxsph_, maxvert_
      real(c_double)  :: xtscor_(maxts_), ytscor_(maxts_), ztscor_(maxts_)
      real(c_double)  :: xsphcor_(maxts_), ysphcor_(maxts_), zsphcor_(maxts_), rsph_(maxts_)
      real(c_double)  :: ar_(maxts_), xe_(maxts_), ye_(maxts_), ze_(maxts_), rin_(maxts_)
      real(c_double)  :: avgArea_, rsolv_, ret_
      integer(c_int)  :: nts_, ntsirr_, nesfp_, addsph_
      logical(c_bool) :: pedra_file_exists
      integer(c_int)  :: pgroup_

      integer(c_int)   :: i
      integer(c_int)   :: error_code
      integer          :: lvpri
      type(point_group) :: pgroup
     
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
      areats = avgArea_
      icesph = 1
      iprpcm = 3
      call get_point_group(pgroup, pgroup_)
      write(lvpri, *) " Point group is ", pgroup%group_name
      rsolv = rsolv_
      omega = 40.0d+00
      fro = 0.7d+00
! ret is the minimum radius of added spheres
      ret = ret_ 
! nesfp is the number of initial spheres.
! nesf is the total number of spheres: initial + added
      nesfp = nesfp_
      do i = 1, nesfp
         xe(i) = xe_(i)
         ye(i) = ye_(i)
         ze(i) = ze_(i)
         rin(i) = rin_(i)
         alpha(i) = 1.0d0
      enddo
      
      nesf = nesfp

      call polyhedra_driver(pgroup, lvpri, error_code)

! Common block dark magic, it will disappear one day...
      nts_ = nts
      ntsirr_ = ntsirr
! Pass the number of added spheres back, to update the GePolCavity
! object in the right way.
      addsph_ = nesf - nesfp
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

      write(lvpri, *) "Error code is ", error_code

      close(lvpri)

      end subroutine generatecavity_cpp
