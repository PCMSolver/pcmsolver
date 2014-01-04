C
C     simple input reader for cavity generator
C     written by Krzysztof Mozgawa, 2010
C
      subroutine generatecavity_cpp(xtscor,ytscor,ztscor,ar,xsphcor,
     &     ysphcor, zsphcor, rsph, nts, nesfp, xe, ye, ze, rin, avgArea,
     $     rsolv, work, lwork) bind(c, name='generatecavity_cpp')

      use, intrinsic :: iso_c_binding    

      implicit none
      real(c_double) :: xtscor(*), ytscor(*), ztscor(*)
      real(c_double) :: xsphcor(*), ysphcor(*), zsphcor(*), rsph(*)
      real(c_double) :: ar(*), xe(*), ye(*), ze(*), rin(*)
      real(c_double) :: avgArea, rsolv, work(*)
      integer(c_int) :: nts, nesfp, lwork

      call makecav_(xtscor,ytscor,ztscor,ar,xsphcor,
     &     ysphcor, zsphcor, rsph, nts, nesfp, xe, ye, ze, rin, avgArea,
     $     rsolv, work, lwork)

      end
