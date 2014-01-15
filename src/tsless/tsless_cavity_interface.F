      subroutine generate_tslesscavity_cpp(xtscor,ytscor,ztscor,ar,
     &     xsphcor, ysphcor, zsphcor, rsph, nts, nesfp, xe, ye, ze, rin,
     &     avgArea, dmin, nord, rsolv, work, lwork) 
     $     bind(c, name='generate_tslesscavity_cpp')

      use, intrinsic :: iso_c_binding    

      implicit none
      real(c_double) :: xtscor(*), ytscor(*), ztscor(*)
      real(c_double) :: xsphcor(*), ysphcor(*), zsphcor(*), rsph(*)
      real(c_double) :: ar(*), xe(*), ye(*), ze(*), rin(*)
      real(c_double) :: avgArea, rsolv, dmin, work(*)
      integer(c_int) :: nts, nesfp, nord, lwork

      call makecav_(xtscor,ytscor,ztscor,ar,xsphcor,
     &     ysphcor, zsphcor, rsph, nts, nesfp, xe, ye, ze, rin, avgArea,
     $     rsolv, work, lwork)

      end subroutine generate_tslesscavity_cpp
